using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace noiSNESs_Brr_Finder
{
    class Program
    {
        private static string _fileName = string.Empty;
        private static string _additionalErrMessage = String.Empty;

        private static int  _minimumSize = 450;
        private static uint _samplerate = 0;

        private static int _minimumStdDeviation   = 80;
        private static int _maximumStdDeviation   = 100;
        private static int _maximumModePercentage = 20;

        /// <summary>
        /// Main.
        /// </summary>
        /// <param name="args">Input arguments.</param>
        static void Main(string[] args)
        {
            if(!parseArgs(args))
            {
                printUsage();
                if (!string.IsNullOrEmpty(_additionalErrMessage))
                    Console.WriteLine(_additionalErrMessage);
                return;
            }

            try
            {
                searchAndExport();
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
            }
        }

        /// <summary>
        /// Using the input data, look for all the pointers to BRRs in the file and
        /// export them in format BRR and WAV, if needed.
        /// </summary>
        private static void searchAndExport()
        {
            byte[] rom = File.ReadAllBytes(_fileName);

            Console.WriteLine("Looking for BRR files...");
            Dictionary<int, int> brrs = Business.findBrrs(rom, _minimumSize, _minimumStdDeviation,
                _maximumStdDeviation, _maximumModePercentage);

            string destinationDirectory =
                Path.GetDirectoryName(_fileName) +
                Path.DirectorySeparatorChar +
                Path.GetFileNameWithoutExtension(_fileName);
            Directory.CreateDirectory(destinationDirectory);
            destinationDirectory += Path.DirectorySeparatorChar;

            int i = 0;
            int maxi = brrs.Count;
            Console.WriteLine("Exporting the " + maxi + " detected BRR files...");
            foreach (var brrFile in brrs)
            {
                printProgressBar(i++, maxi, 80);

                string newName = destinationDirectory + brrFile.Key.ToString("X6") + "_" + brrFile.Value.ToString("X4");
                byte[] brr = rom.Skip(brrFile.Key).Take(brrFile.Value).ToArray();
                File.WriteAllBytes(newName + ".brr", brr);

                if (_samplerate > 0)
                {
                    bool loopedRef = false;
                    byte[] wav = Brr.decodeBRRToWav(brr, ref loopedRef, true, true, 16, 1, _samplerate);
                    File.WriteAllBytes(newName + "_" + _samplerate + ".wav", wav);
                }
            }

            Console.WriteLine();
            Console.WriteLine("Success!");
        }

        #region Usage and arguments

        /// <summary>
        /// Parse the input arguments.
        /// </summary>
        /// <param name="args">Input arguments taken from the command line.</param>
        /// <returns>True if the arguemtns fit the bounds.</returns>
        private static bool parseArgs(string[] args)
        {
            try
            {
                if(args.Length <= 0)
                {
                    _additionalErrMessage = "The file name is a mandatory argument.";
                    return false;
                }

                _fileName = args[0];
                if (!File.Exists(_fileName))
                {
                    _additionalErrMessage = "The file doesn't exist.";
                    return false;
                }
                
                args = args.Skip(1).ToArray();
                while (args.Length > 0)
                {
                    switch (args[0])
                    {
                        case "-m":
                            _minimumSize = int.Parse(args[1]);
                            args = args.Skip(2).ToArray();
                            break;
                        case "-w":
                            _samplerate = uint.Parse(args[1]);
                            args = args.Skip(2).ToArray();
                            break;
                        case "-s":
                            _minimumStdDeviation = int.Parse(args[1]);
                            _maximumStdDeviation = int.Parse(args[2]);
                            _maximumModePercentage = int.Parse(args[3]);
                            args = args.Skip(4).ToArray();
                            break;
                        default:
                            return false;
                    }
                }

                return validateArguments();
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// Validate the input arguments.
        /// </summary>
        /// <returns>True if the arguemtns fit the bounds.</returns>
        private static bool validateArguments()
        {
            if (_minimumSize > 0xFFFF || _minimumSize < 0)
                _additionalErrMessage = "The given value to minimum size is not valid.";

            if (_samplerate < 4000 || _samplerate > 32000)
                _additionalErrMessage = "The samplerate must be a value between 4000 and 32000.";

            if (_minimumStdDeviation < 0 || _minimumStdDeviation > 100)
                _additionalErrMessage = "The minimum standard deviation must be a value between 0 and 100.";

            if (_maximumStdDeviation < 80 || _maximumStdDeviation > 150)
                _additionalErrMessage = "The maximum standard deviation must be a value between 80 and 150.";

            if (_minimumStdDeviation >= _maximumStdDeviation)
                _additionalErrMessage = "The minimum must be lesser than the maximum standard deviation.";

            if (_maximumModePercentage > 100 || _maximumModePercentage < 0)
                _additionalErrMessage = "The maximum mode percentage must be a value between 0 and 100.";

            return string.IsNullOrEmpty(_additionalErrMessage);
        }

        /// <summary>
        /// Print the usage legend in the console
        /// </summary>
        private static void printUsage()
        {
            Console.WriteLine("usage: noiSNESs_Brr_Finder.exe <file_name> [Options]");
            Console.WriteLine("Options:");
            Console.WriteLine("\t-m <size> Minimum size of brr chunk. (Default 450)");
            Console.WriteLine("\t-w <samplerate> Export files as .wav with the given samplerate.");
            Console.WriteLine("\t-s <min_std_deviation max_std_deviation max_mode%> Statistics");
            Console.WriteLine("\t   a brr sample should match to not be considered as a false");
            Console.WriteLine("\t   positives. The default values are 80, 100, and 20.");
            Console.WriteLine();
            Console.WriteLine("Example: noiSNESs_Brr_Finder.exe smw.smc -s 80 100 20 -w 8000");
            Console.WriteLine();
        }

        #endregion

        #region Console helper

        /// <summary>
        /// Print a progress bar in the command line.
        /// </summary>
        /// <param name="progress">Current status of the progress.</param>
        /// <param name="max">Maximum value of the progress.</param>
        /// <param name="barSize">Length of the bar in the command line.</param>
        private static void printProgressBar(int progress, int max, int barSize)
        {
            int i = 0;
            int n = (progress * barSize) / max;

            Console.Write("\r[");
            for (i = 0; i < n; i++)
                Console.Write("=");
            for (; i < barSize; i++)
                Console.Write("·");
            Console.Write("]");
        }

        #endregion
    }
}
