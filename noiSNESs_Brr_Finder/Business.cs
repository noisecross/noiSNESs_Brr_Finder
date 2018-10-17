using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace noiSNESs_Brr_Finder
{
    class Business
    {
        #region BRR detection

        /// <summary>
        /// Find the BRRs contained in a data chunk.
        /// </summary>
        /// <param name="dataChunk">Data chunk to find BRR data in.</param>
        /// <param name="minimumSize">Minimum size to consider a BRR data.</param>
        /// <param name="minimumStdDeviation">Minimum standard deviation to consider valid a BRR data.</param>
        /// <param name="maximumStdDeviation">Maximum standard deviation to consider valid a BRR data.</param>
        /// <param name="maximumModePercentage">Maximum percentage of the mode value to consider valid a BRR data.</param>
        /// <returns>The dictionay of pairs (Address,size) of every found BRR.</returns>
        public static Dictionary<int, int> findBrrs(byte[] dataChunk, int minimumSize = 450,
            int minimumStdDeviation = 80, int maximumStdDeviation = 100, int maximumModePercentage = 20)
        {
            Dictionary<int, int> output = new Dictionary<int, int>();
            List<int> analyzerIndex = new List<int> { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08 };

            while (true)
            {
                //Return when all the indexes are at the end of the file
                if (!analyzerIndex.Any(_ => _ < dataChunk.Length))
                    return output;

                //Analyze the less advanced pointer
                int analyzer = analyzerIndex.IndexOf(analyzerIndex.Min());

                int newIndex = 0;
                if (checkBRR(dataChunk, analyzerIndex[analyzer], out newIndex, minimumSize,
                    minimumStdDeviation, maximumStdDeviation, maximumModePercentage))
                {
                    //Brr found
                    output.Add(analyzerIndex[analyzer], newIndex - analyzerIndex[analyzer]);

                    //Update the indexes which needed it
                    for (int i = 0; i < 9; i++)
                    {
                        int modAnalyzer = (analyzer + i) % 9;

                        if (analyzerIndex[modAnalyzer] < newIndex + i)
                            analyzerIndex[modAnalyzer] = newIndex + i;
                    }
                }
                else
                {
                    analyzerIndex[analyzer] = newIndex;
                }
            }
        }

        /// <summary>
        /// Check if the index is pointing a valid BRR a the data chunk.
        /// </summary>
        /// <param name="dataChunk"></param>
        /// <param name="index">The index to start checking the data.</param>
        /// <param name="newIndex">The new position of the index after performing the analysis.</param>
        /// <param name="minimum_size">Minimum size of a BRR candidate to avoid false positives.</param>
        /// <param name="minimumStdDeviation">Minimum standard deviation to consider valid a BRR data.</param>
        /// <param name="maximumStdDeviation">Maximum standard deviation to consider valid a BRR data.</param>
        /// <param name="maximumModePercentage">Maximum percentage of the mode value to consider valid a BRR data.</param>
        /// <returns>True if the index is pointing a valid BRR in the data chunk.</returns>
        private static bool checkBRR(byte[] dataChunk, int index, out int newIndex, int minimum_size = 450,
            int minimumStdDeviation = 80, int maximumStdDeviation = 100, int maximumModePercentage = 20)
        {
            newIndex = index;

            while (true)
            {
                if (newIndex >= dataChunk.Length)
                    return false;

                if (dataChunk[newIndex] > 0xCF)
                {
                    newIndex += 9;
                    return false;
                }

                if ((dataChunk[newIndex] & 0x01) != 0)
                {
                    //Last byte
                    newIndex += 9;

                    int size = newIndex - index;
                    if (size > minimum_size)
                        return discriminateByModeAndStdDeviation(dataChunk.Skip(index).Take(size).ToList(),
                            minimumStdDeviation, maximumStdDeviation, maximumModePercentage);

                    return false;
                }

                newIndex += 9;
            }
        }

        /// <summary>
        /// Performs an statistical analysis over a byte stream to detect false positives.
        /// </summary>
        /// <param name="brrCandidate">Byte stream candidate to be a valid BRR.</param>
        /// <param name="minimumStdDeviation">Minimum standard deviation to consider valid a BRR data.</param>
        /// <param name="maximumStdDeviation">Maximum standard deviation to consider valid a BRR data.</param>
        /// <param name="maximumModePercentage">Maximum percentage of the mode value to consider valid a BRR data.</param>
        /// <returns>True if the input stream is a valid BRR.</returns>
        private static bool discriminateByModeAndStdDeviation(List<byte> brrCandidate,
            int minimumStdDeviation = 80, int maximumStdDeviation = 100, int maximumModePercentage = 20)
        {
            try
            {
                byte mode = brrCandidate.GroupBy(item => item).OrderByDescending(item => item.Count()).First().Key;
                double modePercentage = 100 * brrCandidate.Count(_ => _ == mode) / (double)brrCandidate.Count;
                if (modePercentage > maximumModePercentage)
                    return false;

                double average = brrCandidate.Average(_ => _);
                double stdDeviation = Math.Sqrt(brrCandidate.Average(_ => Math.Pow(_ - average, 2)));

                return (stdDeviation > minimumStdDeviation) && (stdDeviation < maximumStdDeviation);
            }
            catch
            {
                return false;
            }
        }

        #endregion
    }
}
