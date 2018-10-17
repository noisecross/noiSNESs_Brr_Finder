using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace noiSNESs_Brr_Finder
{
    public sealed class Brr
    {
        //------------------------------------------
        // IMPORTANT NOTE
        // This class its a copy and translation to
        // C# made from the code written by 
        // Bregalad in his BRRtools. All the credit
        // from this class goes to him.
        //------------------------------------------

        public static readonly int BRR_CHUNK_SIZE = 9;

        public static byte[] decodeBRR(byte[] brr, ref bool ptrLooped)
        {
            List<Int16> rawSamples = decodeBRRToInt16(brr, ref ptrLooped);
            return rawSamples.SelectMany(BitConverter.GetBytes).ToArray();
        }

        public static void playBrr(byte[] brr, uint sampleRate)
        {
            bool loopedRef = false;

            byte[] wav = decodeBRRToWav(brr, ref loopedRef, true, true, 16, 1, sampleRate);
            using (MemoryStream msw = new MemoryStream(wav))
            {
                // Construct the sound player
                System.Media.SoundPlayer playerW = new System.Media.SoundPlayer(msw);
                playerW.Play();
            }
        }

        //This process is specific to play STAR OCEAN samples.
        public static void playSOStreamedBrr(byte[] brr, uint sampleRate)
        {
            bool loopedRef = false;

            byte[] wav = decodeSOStreamBRRToWav(ref brr, ref loopedRef, true, true, 16, 1, sampleRate);
            using (MemoryStream msw = new MemoryStream(wav))
            {
                // Construct the sound player
                System.Media.SoundPlayer playerW = new System.Media.SoundPlayer(msw);
                playerW.Play();
            }
        }

        private static void gaussFilter(ref List<Int16> buffer)
        {
            int length = buffer.Count;

            int prev = (372 + 1304) * buffer[0] + 372 * buffer[1]; // First sample
            for (int i = 1; i < length - 1; ++i)
            {
                int k0 = 372 * (buffer[i - 1] + buffer[i + 1]);
                int k = 1304 * buffer[i];
                buffer[i - 1] = Convert.ToInt16(prev / 2048);
                prev = k0 + k;
            }
            int last = 372 * buffer[length - 2] + (1304 + 372) * buffer[length - 1];
            buffer[length - 2] = Convert.ToInt16(prev / 2048);
            buffer[length - 1] = Convert.ToInt16(last / 2048);
        }

        #region Export brr

        public static List<Int16> decodeBRRToInt16(byte[] brr, ref bool ptrLooped, bool checkMultiplicy=true)
        {
            List<Int16> rawSamples = new List<Int16>();

            //Check tamaño es múltiplo de 9
            if (checkMultiplicy && (brr.Length % 9) != 0)
                return rawSamples;

            Int32[] prev = { 0, 0 };
            int decodedSize = 0;

            while (decodedSize + BRR_CHUNK_SIZE <= brr.Length)
            {
                byte flags = brr[decodedSize];
                bool chunkEnd = (flags & 1) != 0;
                bool chunkLoop = (flags & 2) != 0;
                byte filter = Convert.ToByte((flags >> 2) & 3);
                byte range = Convert.ToByte(flags >> 4);
                bool validRange = (range <= 0x0c);

                Int32 S1 = prev[0];
                Int32 S2 = prev[1];

                for (int byteIndex = 0; byteIndex < 8; byteIndex++)
                {
                    sbyte sample1 = (sbyte)brr[decodedSize + 1 + byteIndex];
                    sbyte sample2 = (sbyte)(sample1 << 4);
                    sample1 >>= 4;
                    sample2 >>= 4;

                    for (int nybble = 0; nybble < 2; nybble++)
                    {
                        Int32 output;

                        output = (nybble != 0) ? (Int32)sample2 : (Int32)sample1;
                        output = validRange ? ((output << range) >> 1) : (output & ~0x7FF);

                        switch (filter)
                        {
                            case 0: // Direct
                                break;

                            case 1: // 15/16
                                output += S1 + ((-S1) >> 4);
                                break;

                            case 2: // 61/32 - 15/16
                                output += (S1 << 1) + ((-((S1 << 1) + S1)) >> 5) - S2 + (S2 >> 4);
                                break;

                            case 3: // 115/64 - 13/16
                                output += (S1 << 1) + ((-(S1 + (S1 << 2) + (S1 << 3))) >> 6) - S2 + (((S2 << 1) + S2) >> 4);
                                break;
                        }

                        output = sclip15(sclamp16(output));

                        S2 = S1;
                        S1 = output;

                        rawSamples.Add(Convert.ToInt16(output << 1));
                    }
                }

                prev[0] = S1;
                prev[1] = S2;

                decodedSize += BRR_CHUNK_SIZE;

                if (chunkEnd)
                {
                    ptrLooped = chunkLoop;
                    break;
                }
            }

            return rawSamples;
        }

        public static byte[] decodeBRRToWav(byte[] brr, ref bool ptrLooped, bool addWavHeader = true, bool applyGaussFilter = true, Int16 bitwidth = 16, Int16 channels = 1, UInt32 samplerate = 32000)
        {
            List<Int16> rawSamples = decodeBRRToInt16(brr, ref ptrLooped);
            if (rawSamples.Count <= 0)
                return new byte[] { };

            if (applyGaussFilter)
                gaussFilter(ref rawSamples);

            if (!addWavHeader)
                return rawSamples.SelectMany(BitConverter.GetBytes).ToArray();

            byte[] wavBytes = Brr.getWavHeader(rawSamples.Count, 16, 1, samplerate).ToArray();
            wavBytes = wavBytes.Concat(rawSamples.SelectMany(BitConverter.GetBytes).ToArray()).ToArray();
            return wavBytes;
        }

        //This process is specific to decode STAR OCEAN samples.
        public static byte[] decodeSOStreamBRRToWav(ref byte[] brr, ref bool ptrLooped, bool addWavHeader = true, bool applyGaussFilter = true, Int16 bitwidth = 16, Int16 channels = 1, UInt32 samplerate = 32000)
        {
            List<Int16> rawSamples = new List<Int16>();

            int nOfBrrBytes = brr[7] + (brr[8] & 0x3F) * 0x0100;
            brr = brr.Skip(9).ToArray();
            brr = brr.Take(nOfBrrBytes).ToArray();
            nOfBrrBytes = (nOfBrrBytes / 9) * 16;

            while (rawSamples.Count < nOfBrrBytes)
            {
                int startIndex = (rawSamples.Count * 9) / 16;
                rawSamples.AddRange(decodeBRRToInt16(brr.Skip(startIndex).ToArray(), ref ptrLooped, false));
            }
            
            if (rawSamples.Count <= 0)
                return new byte[] { };

            if (applyGaussFilter)
                gaussFilter(ref rawSamples);

            if (!addWavHeader)
                return rawSamples.SelectMany(BitConverter.GetBytes).ToArray();

            byte[] wavBytes = Brr.getWavHeader(rawSamples.Count, 16, 1, samplerate).ToArray();
            wavBytes = wavBytes.Concat(rawSamples.SelectMany(BitConverter.GetBytes).ToArray()).ToArray();
            return wavBytes;
        }

        public static List<byte> getWavHeader(Int32 nSamples, Int16 bitwidth = 16, Int16 channels = 1, UInt32 samplerate = 32000)
        {
            List<byte> header = new List<byte>();
            Int16 bytesPerSample = Convert.ToInt16(bitwidth / 8);
            UInt32 fileSize = Convert.ToUInt32(36 + (bytesPerSample * nSamples));

            //Header      36 bytes
            header.AddRange(Encoding.ASCII.GetBytes("RIFF"));                                        //00
            header.AddRange(BitConverter.GetBytes(fileSize));                                        //04
            header.AddRange(Encoding.ASCII.GetBytes("WAVE"));                                        //08
            header.AddRange(Encoding.ASCII.GetBytes("fmt "));                                        //12
            header.AddRange(BitConverter.GetBytes((Int32)16));                                       //16
            header.AddRange(BitConverter.GetBytes((Int16)1));                                        //18
            header.AddRange(BitConverter.GetBytes(channels));                                        //20
            header.AddRange(BitConverter.GetBytes(samplerate));                                      //22
            header.AddRange(BitConverter.GetBytes((Int32)(samplerate * bytesPerSample * channels))); //26
            header.AddRange(BitConverter.GetBytes((Int16)(bytesPerSample * channels)));              //30
            header.AddRange(BitConverter.GetBytes(bitwidth));                                        //32

            //Data header  8 bytes
            header.AddRange(Encoding.ASCII.GetBytes("data"));                  //36
            header.AddRange(BitConverter.GetBytes(nSamples * bytesPerSample)); //40
            //44
            return header;
        }

        #endregion

        #region Import Wav as Brr

        public static byte[] encodeBrr(WavFileHeader hdr, Int16[] samples, BrrEncoderOptions options)
        {
            List<byte> output = new List<byte>();

            byte loopFlag = 0x00;
            UInt16 loopStart = 0;
            uint newLoopsize = 0;

            if (options.loopAddress >= 0)
            {
                loopFlag = 0x02; // 0x02 if loop flag is active.
                loopStart = Convert.ToUInt16(options.loopAddress);
            }

            // Optional truncation of input sample
            UInt32 samplesLength = hdr.dataSize / hdr.blockAlign; // FIXME comprobar qué es esto ¿no se puede utilizar samples.Length?            

            // Adjust amplitude in function of amount of channels
            samples = adjustAmplitudeOfWavSamples(samples, samplesLength, hdr.chans, hdr.bitsPerSample, options);

            // Apply resampling if needed
            applyResampling(ref samples, ref samplesLength, ref newLoopsize, hdr.sampleRate, loopStart, options);

            // Apply trebble boost filter (gussian lowpass compensation) if requested by user
            if (options.trebleBoost)
                samples = trebleBoostFilter(samples, samplesLength);

            //Add zeroes at the beginning if the Amount of PCM samples isn't a multiple of 16
            int zeroesToAdd = samples.Length % 16;
            if (zeroesToAdd != 0)
                samples = Enumerable.Repeat((Int16)0, zeroesToAdd).Concat(samples).ToArray();

            //Begin the actual BRR encoding

            //Initialization needed if any of the first 16 samples isn't zero
            bool needsInitialBrrBlock = false;
            if (samples.Length >= 16)
            {
                for (int i = 0; i < 16 && !needsInitialBrrBlock; ++i)
                    needsInitialBrrBlock |= samples[i] != 0;
            }
            else
                needsInitialBrrBlock = true;

            if (needsInitialBrrBlock)
                output = new List<byte>() { loopFlag, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

            //Encode BRR

            //Set the truncate size
            if (options.truncateLen != 0 && output.Count > 0)
                options.truncateLen -= 16;

            if (options.truncateLen != 0 && (options.truncateLen != samplesLength))
                samplesLength = options.truncateLen;

            Int16 p1 = 0;
            Int16 p2 = 0;
            BrrEncoderStatistics statistics = new BrrEncoderStatistics();
            for (int n = 0; n < samplesLength; n += 16)
            {
                byte[] brrChunk = new byte[] { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

                //Encode BRR block, tell the encoder if we're at loop point (if loop is enabled), and if we're at end point
                bool isLoopPoint = options.fixLoopEn && (n == (samplesLength - newLoopsize));
                bool isEndPoint = n == samplesLength - 16;

                Int16[] samplesToCheck = samples.Skip(n).ToArray();
                if (samplesToCheck.Length < 16)
                    samplesToCheck = samplesToCheck.Concat(new Int16[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }).ToArray();

                ADPCMBlockMash(samples.Skip(n).ToArray(), isLoopPoint, isEndPoint, options, statistics, ref p1, ref p2, ref brrChunk);

                //Set the loop flag if needed
                brrChunk[0] |= loopFlag;
                output.AddRange(brrChunk);
            }

            //HACKME recover this
            //    if(fix_loop_en)
            //    {
            //        unsigned int k = samples_length - (initial_block ? new_loopsize - 16 : new_loopsize);
            //        printf("Position of the loop within the BRR sample : %u samples = %u BRR blocks.\n", k, k/16);
            //    }

            //    for(int i=0; i<4; i++)
            //        if (FIRstats[i]>0) printf("Filter %u used on %u blocks.\n", i, FIRstats[i]);

            return output.ToArray();
        }

        public static WavFileHeader openWavFile(string fileName, out Int16[] wavSampleData)
        {
            WavFileHeader hdr = new WavFileHeader();
            wavSampleData = new Int16[] { };
            byte[] wavData;

            try
            {
                wavData = File.ReadAllBytes(fileName);
                if (wavData.Length < 36)
                {
                    hdr.invalidHeaderMessage = "Error : Input file in incompatible format.";
                    return hdr;
                }
            }
            catch (Exception ex)
            {
                hdr.invalidHeaderMessage = "Error : " + ex.Message;
                return hdr;
            }

            return openWavFile(wavData, out wavSampleData);
        }

        public static WavFileHeader openWavFile(byte[] inputWavData, out Int16[] wavSampleData)
        {
            WavFileHeader hdr = new WavFileHeader();
            List<byte> wavData = inputWavData.ToList();
            wavSampleData = new Int16[] { };

            hdr.chunkID = wavData.GetRange(0, 4).ToArray();                                // Should be 'RIFF'
            hdr.chunkSize = BitConverter.ToUInt32(wavData.GetRange(4, 4).ToArray(), 0);
            hdr.waveStr = wavData.GetRange(8, 4).ToArray();                                // Should be 'WAVE'
            hdr.sc1Id = wavData.GetRange(12, 4).ToArray();                                 // Should be 'fmt '
            hdr.sc1size = BitConverter.ToUInt32(wavData.GetRange(16, 4).ToArray(), 0);     // Should be at least 16
            hdr.audioFormat = BitConverter.ToUInt16(wavData.GetRange(20, 2).ToArray(), 0); // Should be 1 for PCM
            hdr.chans = BitConverter.ToUInt16(wavData.GetRange(22, 2).ToArray(), 0);       // 1 for mono, 2 for stereo, etc
            hdr.sampleRate = BitConverter.ToUInt32(wavData.GetRange(24, 4).ToArray(), 0);
            hdr.byteRate = BitConverter.ToUInt32(wavData.GetRange(28, 4).ToArray(), 0);
            hdr.blockAlign = BitConverter.ToUInt16(wavData.GetRange(32, 2).ToArray(), 0);
            hdr.bitsPerSample = BitConverter.ToUInt16(wavData.GetRange(34, 2).ToArray(), 0);

            //Validate header
            string chunkID = Encoding.UTF8.GetString(hdr.chunkID, 0, 4);
            string waveStr = Encoding.UTF8.GetString(hdr.waveStr, 0, 4);
            string sc1Id = Encoding.UTF8.GetString(hdr.sc1Id, 0, 4);

            if (chunkID != "RIFF")
                hdr.invalidHeaderMessage = "Error : Input file in unsupported format : \"RIFF\" block missing.";
            else if (waveStr != "WAVE")
                hdr.invalidHeaderMessage = "Error : Input file in unsupported format : \"WAVE\" block missing.";
            else if (sc1Id != "fmt ")
                hdr.invalidHeaderMessage = "Error : Input file in unsupported format : \"fmt \" block missing.";
            else if (hdr.sc1size < 0x10 || hdr.audioFormat != 1)
                hdr.invalidHeaderMessage = "Input file in unsupported format : file must be uncompressed PCM.";
            else if (hdr.byteRate != hdr.sampleRate * hdr.chans * hdr.bitsPerSample / 8)
                hdr.invalidHeaderMessage = "Byte rate in input file is set incorrectly.";
            else if (hdr.blockAlign != hdr.bitsPerSample * hdr.chans / 8)
                hdr.invalidHeaderMessage = "Block align in input file is set incorrectly.";
            else if (hdr.bitsPerSample != 8 && hdr.bitsPerSample != 16)
                hdr.invalidHeaderMessage = "Error : unsupported amount of bits per sample (8 or 16 are supported)";
            else
                hdr.validHeader = true;

            if (!hdr.validHeader)
                return hdr;

            if (hdr.chans != 1)
                hdr.invalidHeaderMessage = "Input is multi-channel : Will automatically be converted to mono.";

            //Validate data header
            UInt32 index = 0x24 + hdr.sc1size - 0x10;
            while (true)
            {
                try
                {
                    string subHeaderName = Encoding.UTF8.GetString(wavData.GetRange((int)index, 4).ToArray(), 0, 4);
                    UInt32 subHeaderSize = BitConverter.ToUInt32(wavData.GetRange((int)index + 4, 4).ToArray(), 0);
                    index += 8;

                    if (subHeaderName == "data")
                    {
                        hdr.dataIndex = index;
                        hdr.dataSize = subHeaderSize;

                        List<Int16> buffer = new List<Int16>();
                        using (MemoryStream ms = new MemoryStream(wavData.GetRange((int)index, (int)subHeaderSize).ToArray()))
                        using (BinaryReader br = new BinaryReader(ms))
                        {
                            try
                            {
                                while (true)
                                    buffer.Add(br.ReadInt16());
                            }
                            catch (EndOfStreamException) { }
                        }

                        wavSampleData = buffer.ToArray();
                        break;
                    }

                    index += subHeaderSize;
                }
                catch
                {
                    hdr.validHeader = false;
                    hdr.invalidHeaderMessage = "End of file reached without finding a correct \"data\" chunk.";
                    return hdr;
                }
            }

            //Everything is right. Return hdr.            
            return hdr;
        }

        #endregion

        #region private functions

        private static int getBrrPrediction(byte filter, Int16 p1, Int16 p2)
        {
            int p;

            //Different formulas for 4 filters
            switch (filter)
            {
                case 0:
                    return 0;

                case 1:
                    p = p1;
                    p -= p1 >> 4;
                    return p;

                case 2:
                    p = p1 << 1;
                    p += (-(p1 + (p1 << 1))) >> 5;
                    p -= p2;
                    p += p2 >> 4;
                    return p;

                case 3:
                    p = p1 << 1;
                    p += (-(p1 + (p1 << 2) + (p1 << 3))) >> 6;
                    p -= p2;
                    p += (p2 + (p2 << 1)) >> 4;
                    return p;

                default:
                    return 0;
            }
        }

        // Convert a block from PCM to BRR
        // Returns the squared error between original data and encoded data
        // If "is_end_point" is true, the predictions p1/p2 at loop are also used in caluclating the error (depending on filter at loop)
        private static double ADPCMMash(int shiftamount, byte filter, Int16[] PCMData, bool write, bool isEndPoint, BrrEncoderOptions options, BrrEncoderStatistics statistics, ref Int16 p1, ref Int16 p2, ref byte[] brr)
        {
            double d2 = 0.0;
            Int16 l1 = p1;
            Int16 l2 = p2;
            int step = 1 << shiftamount;

            int vlin, d, da, dp, c;

            if (PCMData.Length < 16)
                return 0;

            for (int i = 0; i < 16; ++i)
            {
                //Make linear prediction for next sample
                //vlin = (v0 * iCoef[0] + v1 * iCoef[1]) >> 8;
                vlin = getBrrPrediction(filter, l1, l2) >> 1;
                d = (PCMData[i] >> 1) - vlin; //difference between linear prediction and current sample
                da = Math.Abs(d);

                // Take advantage of wrapping
                if (options.wrapEn && da > 16384 && da < 32768)
                {
                    d = d - 32768 * (d >> 24);
                    //if(write) printf("Caution : Wrapping was used.\n");
                }

                dp = d + (step << 2) + (step >> 2);
                c = 0;
                if (dp > 0)
                {
                    if (step > 1)
                        c = dp / (step / 2);
                    else
                        c = dp * 2;
                    if (c > 15)
                        c = 15;
                }
                c -= 8;
                dp = (c << shiftamount) >> 1; // quantized estimate of samp - vlin
                                              // edge case, if caller even wants to use it
                if (shiftamount > 12)
                    dp = (dp >> 14) & ~0x7FF;
                c &= 0x0f; // mask to 4 bits

                l2 = l1; //shift history */
                l1 = (Int16)(sclamp16(vlin + dp) * 2);

                d = PCMData[i] - l1;
                d2 += (double)d * d; // update square-error

                if (write) // If we want output, put it in proper place
                {
                    brr[1 + (i >> 1)] |= ((i & 1) != 0) ? (byte)c : (byte)(c << 4);
                }
            }

            // Also account for history points when looping is enabled & filters used
            if (!isEndPoint)
                d2 /= 16.0;
            else
            {
                switch (statistics.filterAtLoop)
                {
                    case 0:
                        d2 /= 16.0;
                        break;

                    case 1: // Filter 1
                        d = l1 - statistics.p1AtLoop;
                        d2 += (double)d * d;
                        d2 /= 17.0;
                        break;

                    default: // Filters 2 & 3
                        d = l1 - statistics.p1AtLoop;
                        d2 += (double)d * d;
                        d = l2 - statistics.p2AtLoop;
                        d2 += (double)d * d;
                        d2 /= 18.0;
                        break;
                }
            }

            // when generating real output, we want to return these
            if (write)
            {
                p1 = l1;
                p2 = l2;

                //Set the end bit if we're on the last block
                brr[0] = Convert.ToByte((shiftamount << 4) | (filter << 2));
                if (isEndPoint)
                    brr[0] |= 1;
            }
            return d2;
        }

        // Encode a ADPCM block using brute force over filters and shift amounts
        private static void ADPCMBlockMash(Int16[] PCMData, bool isLoopPoint, bool isEndPoint, BrrEncoderOptions options, BrrEncoderStatistics statistics, ref Int16 p1, ref Int16 p2, ref byte[] brr)
        {
            byte smin = byte.MaxValue;
            byte kmin = byte.MaxValue;
            double dmin = double.MaxValue;

            for (byte s = 0; s < 13; ++s)
            {
                for (byte k = 0; k < 4; ++k)
                {
                    if (!options.FIRen[k])
                        continue;
                    double d = ADPCMMash(s, k, PCMData, false, isEndPoint, options, statistics, ref p1, ref p2, ref brr);
                    if (d >= dmin)
                        continue;

                    kmin = k;		//Memorize the filter, shift values with smaller error
                    dmin = d;
                    smin = s;
                }
            }

            if (isLoopPoint)
            {
                statistics.filterAtLoop = kmin;
                statistics.p1AtLoop = p1;
                statistics.p2AtLoop = p2;
            }

            ADPCMMash(smin, kmin, PCMData, true, isEndPoint, options, statistics, ref p1, ref p2, ref brr);
            statistics.FIRstats[kmin]++;
        }

        // This function applies a treble boosting filter that compensates the gauss lowpass filter
        private static Int16[] trebleBoostFilter(Int16[] samples, uint length)
        {
            // Tepples' coefficient multiplied by 0.6 to avoid overflow in most cases
            double[] coefs = new double[] { 0.912962, -0.16199, -0.0153283, 0.0426783, -0.0372004, 0.023436, -0.0105816, 0.00250474 };

            Int16[] output = new Int16[length];

            for (int i = 0; i < length; ++i)
            {
                double acc = samples[i] * coefs[0];

                for (int k = 7; k > 0; --k)
                {
                    acc += coefs[k] * ((i + k < length) ? samples[i + k] : samples[length - 1]);
                    acc += coefs[k] * ((i - k >= 0) ? samples[i - k] : samples[0]);
                }

                output[i] = (Int16)acc;
            }

            return output;
        }

        private static Int16[] adjustAmplitudeOfWavSamples(Int16[] samples, UInt32 samplesLength, int chans, int bitsPerSample, BrrEncoderOptions options)
        {
            List<Int16> output = new List<Int16>();

            // Adjust amplitude in function of amount of channels
            options.amplAdjust /= chans;

            int sample;
            if (bitsPerSample == 8)
            {
                List<byte> input8 = samples.SelectMany(BitConverter.GetBytes).ToList();

                for (int i = 0; i < samplesLength; ++i)
                {
                    //fread(in8_chns, 1, hdr.chans, inwav);   // Read single sample on CHANS channels at a time
                    List<byte> in8Chns = input8.GetRange(i * chans, chans);
                    sample = 0;
                    for (int ch = 0; ch < chans; ++ch)      // Average samples of all channels
                        sample += in8Chns[ch] - 0x80;
                    output.Add((Int16)((sample << 8) * options.amplAdjust));
                }
            }
            else if (bitsPerSample == 16)
            {
                List<Int16> input16 = samples.ToList();

                for (int i = 0; i < samplesLength; ++i)
                {
                    //fread(in16_chns, 2, hdr.chans, inwav);
                    List<Int16> in16Chns = input16.GetRange(i * chans, chans);
                    sample = 0;

                    for (int ch = 0; ch < chans; ++ch)
                        sample += in16Chns[ch];
                    output.Add((Int16)(sample * options.amplAdjust));
                }
            }
            //else Error : unsupported amount of bits per sample (8 or 16 are supported)

            return output.ToArray();
        }

        private static void applyResampling(ref Int16[] samples, ref uint samplesLength, ref uint newLoopsize, uint sourceSamplerate, uint loopStart, BrrEncoderOptions options)
        {
            uint targetLength;

            // Set resample factor if auto samplerate mode
            if (options.targetSamplerate != 0)
                targetLength = (samplesLength * options.targetSamplerate) / sourceSamplerate;
            else
                targetLength = Convert.ToUInt32(samplesLength / options.ratio);

            if (options.fixLoopEn)
            {
                uint loopsize = ((samplesLength - loopStart) * targetLength) / samplesLength;

                // Adjust resampling
                newLoopsize = ((loopsize + 15) / 16) * 16; // New loopsize is the multiple of 16 that comes after loopsize
                targetLength = (targetLength * newLoopsize) / loopsize;
            }

            samples = resample(samples, samplesLength, targetLength, options);
            samplesLength = targetLength;
        }

        private static Int16[] resample(Int16[] samples, uint samplesLength, uint outLength, BrrEncoderOptions options)
        {
            Int16[] output;

            char type = options.resampleType;
            options.ratio = (double)samplesLength / (double)outLength;

            switch (type)
            {
                case 'l': // Linear interpolation
                    output = WaveformsInterpolation.linearInterpolation(samples, options.ratio, outLength);
                    break;
                case 's': // Sine interpolation
                    output = WaveformsInterpolation.sineInterpolation(samples, options.ratio, outLength, samplesLength);
                    break;
                case 'c': // Cubic interpolation
                    output = WaveformsInterpolation.cubicInterpolation(samples, options.ratio, outLength, samplesLength);
                    break;
                case 'b': // Bandlimited interpolation
                    output = WaveformsInterpolation.bandLimitedInterpolation(samples, options.ratio, outLength, samplesLength);
                    break;
                case 'n': // No interpolation
                default:  // Undefined
                    output = new Int16[outLength];
                    for (int i = 0; i < outLength; ++i)
                        output[i] = samples[Convert.ToInt32(Math.Min(Math.Floor(i * options.ratio), samplesLength - 1))];
                    break;
            };

            return output.ToArray();
        }

        private static Int32 sclip15(Int32 x)
        {
            if ((x & 16384) != 0)
                return x | ~16383;
            else
                return x & 16383;
        }

        private static Int32 sclamp16(Int32 x)
        {
            if (x > 32767)
                return 32767;

            return (x < -32768) ? -32768 : x;
        }

        private static class WaveformsInterpolation
        {
            public static Int16[] linearInterpolation(Int16[] samples, double ratio, uint outLength)
            {
                //Linear interpolation
                Int16[] output = new Int16[outLength];

                int a = 0;
                for (int i = 0; i < outLength; ++i)
                {
                    a = (int)(i * ratio);     //Whole part of index
                    double b = i * ratio - a; //Fractional part of index
                    if ((a + 1) >= samples.Length)
                        output[(int)outLength - 1] = samples[a];
                    else
                        output[i] = Convert.ToInt16((1 - b) * samples[a] + b * samples[a + 1]);
                }


                return output;
            }

            public static Int16[] sineInterpolation(Int16[] samples, double ratio, uint outLength, uint samplesLength)
            {
                Int16[] output = new Int16[outLength];

                //Sine interpolation
                for (int i = 0; i < outLength; ++i)
                {
                    int a = (int)(i * ratio);
                    double b = i * ratio - a;
                    double c = (1.0 - Math.Cos(b * Math.PI)) / 2.0;

                    if ((a + 1) == samplesLength)
                        output[i] = samples[a]; //This used only for the last sample
                    else
                        output[i] = Convert.ToInt16((1 - c) * samples[a] + c * samples[a + 1]);
                }

                return output;
            }

            public static Int16[] cubicInterpolation(Int16[] samples, double ratio, uint outLength, uint samplesLength)
            {
                Int16[] output = new Int16[outLength];

                //Cubic interpolation
                for (int i = 0; i < outLength; ++i)
                {
                    Int32 a = (Int32)Math.Min(i * ratio, samplesLength - 1);

                    Int16 s0 = (a == 0) ? samples[0] : samples[a - 1];
                    Int16 s1 = samples[a];
                    Int16 s2 = (a + 1 >= samplesLength) ? samples[samplesLength - 1] : samples[a + 1];
                    Int16 s3 = (a + 2 >= samplesLength) ? samples[samplesLength - 1] : samples[a + 2];

                    double a0 = s3 - s2 - s0 + s1;
                    double a1 = s0 - s1 - a0;
                    double a2 = s2 - s0;
                    double b = i * ratio - a;
                    double b2 = b * b;
                    double b3 = b2 * b;
                    output[i] = Convert.ToInt16(b3 * a0 + b2 * a1 + b * a2 + s1);
                }

                return output;
            }

            public static Int16[] bandLimitedInterpolation(Int16[] samples, double ratio, uint outLength, uint samplesLength)
            {
                const int FIR_ORDER = 15;
                Int16[] output = new Int16[outLength];

                // Bandlimited interpolation

                // Antialisaing pre-filtering
                if (ratio > 1.0)
                {
                    Int16[] samplesAntialiased = new Int16[samplesLength];
                    double[] firCoefs = new double[FIR_ORDER + 1];

                    // Compute FIR coefficients
                    for (int k = 0; k <= FIR_ORDER; ++k)
                        firCoefs[k] = sinc(k / ratio) / ratio;

                    // Apply FIR filter to samples
                    for (int i = 0; i < samplesLength; ++i)
                    {
                        double acc = samples[i] * firCoefs[0];
                        for (int k = FIR_ORDER; k > 0; --k)
                        {
                            acc += firCoefs[k] * ((i + k < samplesLength) ? samples[i + k] : samples[samplesLength - 1]);
                            acc += firCoefs[k] * ((i - k >= 0) ? samples[i - k] : samples[0]);
                        }
                        samplesAntialiased[i] = (Int16)acc;
                    }

                    samples = samplesAntialiased;
                }

                // Actual resampling using sinc interpolation
                for (int i = 0; i < outLength; ++i)
                {
                    double a = i * ratio;
                    double acc = 0.0;
                    for (int j = (int)a - FIR_ORDER; j <= (int)a + FIR_ORDER; ++j)
                    {
                        Int16 sample;
                        if (j >= 0)
                            if (j < samplesLength)
                                sample = samples[j];
                            else
                                sample = samples[samplesLength - 1];
                        else
                            sample = samples[0];

                        acc += sample * sinc(a - j);
                    }
                    output[i] = (Int16)acc;
                }

                return output;
            }

            private static double sinc(double x)
            {
                if (x == 0.0)
                    return 1.0;
                else
                    return Math.Sin(Math.PI * x) / (Math.PI * x);
            }
        }

        #endregion
    }

    public sealed class BrrEncoderOptions
    {
        public double amplAdjust = 1.0;     // Adjusting amplitude
        public double ratio = 1.0;          // Resampling factor (range (0..4]) 
        public Int32 loopAddress = -1;      //
        public UInt16 targetSamplerate = 0; // Output sample rate (0 = don't change)
        public bool fixLoopEn = false;      // True if fixed loop is activated
        public UInt16 truncateLen = 0;      // Point at which input wave will be truncated (if = 0, input wave is not truncated)
        public bool trebleBoost = false;    // True to boost the treble in the output BRR
        public bool wrapEn = true;          // Enable wrapping (encoded sample would not be compatible with old SPC players)
        public char resampleType = 'l';     // Resampling type (n = nearest neighboor, l = linear, c = cubic, s = sine, b = bandlimited)

        public bool[] FIRen = new bool[4] { true, true, true, true }; // Which BRR filters are enabled
    }

    public sealed class BrrEncoderStatistics
    {
        public Int16 p1AtLoop = 0;
        public Int16 p2AtLoop = 0;
        public byte filterAtLoop = 0;
        public UInt16[] FIRstats = new UInt16[4] { 0, 0, 0, 0 }; // Statistics on BRR filter usage
    }

    public sealed class WavFileHeader
    {
        public byte[] chunkID;     // Should be 'RIFF'
        public UInt32 chunkSize;
        public byte[] waveStr;     // Should be 'WAVE'
        public byte[] sc1Id;       // Should be 'fmt '
        public UInt32 sc1size;     // Should be at least 16
        public UInt16 audioFormat; // Should be 1 for PCM
        public UInt16 chans;       // 1 for mono, 2 for stereo, etc
        public UInt32 sampleRate;
        public UInt32 byteRate;
        public UInt16 blockAlign;
        public UInt16 bitsPerSample;

        public UInt32 dataIndex;
        public UInt32 dataSize;

        public bool validHeader = false;
        public string invalidHeaderMessage = string.Empty;
    }
}
