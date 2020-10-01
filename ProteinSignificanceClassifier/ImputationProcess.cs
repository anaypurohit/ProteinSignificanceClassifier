using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Text;

namespace ProteinSignificanceClassifier
{
    /// <summary>
    /// Class used for imputing missing intensity values for each protein
    /// </summary>
    public class ImputationProcess
    {
        /// <summary>
        /// Detects protein's missing intensity values and increase the number of such missing values found for each sample
        /// </summary>
        private void CalculateMissingIntensityValuesInSample(ProteinRowInfo proteinRowInfo, double[] missingFactor,
            List<string> samplesFileNames)
        {
            Dictionary<string, double> samplesintensityData = proteinRowInfo.getSamplesIntensityValues();
            for (int i = 0; i < samplesFileNames.Count; i++)
            {
                // detected protein's sample with missing intensity value
                if (samplesintensityData[samplesFileNames[i]] == 0)
                {
                    // increased missing intensity values counter for the sample
                    missingFactor[i]++;
                }
            }
        }

        /// <summary>
        /// Uses protein intensity information to add and eventually compute the sum of all intensity values for each sample. 
        /// Similarly calculates total number of valid intensity values for each sample
        /// </summary>
        private void CalculateNumberAndSumOfIntensityValuesOfSample(ProteinRowInfo proteinRowInfo, double[] sampleAllIntensityValuesSum,
            List<string> samplesFileNames, int[] numberOfIntensityValuesInSample)
        {
            Dictionary<string, double> samplesintensityData = proteinRowInfo.getSamplesIntensityValues();
            for (int i = 0; i < samplesFileNames.Count; i++)
            {
                if (samplesintensityData[samplesFileNames[i]] == 0)
                {
                    continue;
                }

                // increase the total intensity sum for the sample
                sampleAllIntensityValuesSum[i] = sampleAllIntensityValuesSum[i] + samplesintensityData[samplesFileNames[i]];

                // increase number of valid intensity values for the sample
                numberOfIntensityValuesInSample[i]++;
            }
        }

        /// <summary>
        /// Uses protein intensity information to compute the numerator of standard deviation of all intensity values for each sample. 
        /// </summary>
        private void CalculateSampleStandardDeviationNumerator(ProteinRowInfo proteinRowInfo, double[] samplesStandardDeviationNumerators,
            List<string> samplesFileNames, double[] samplesMeanIntensityValue)
        {
            Dictionary<string, double> samplesintensityData = proteinRowInfo.getSamplesIntensityValues();
            for (int i = 0; i < samplesFileNames.Count; i++)
            {
                if (samplesintensityData[samplesFileNames[i]] == 0) continue;

                // increases the numerator component for computing standard deviation of intensity values of sample
                samplesStandardDeviationNumerators[i] = samplesStandardDeviationNumerators[i] +
                    (samplesintensityData[samplesFileNames[i]] - samplesMeanIntensityValue[i]) * (samplesintensityData[samplesFileNames[i]] - samplesMeanIntensityValue[i]);
            }
        }

        /// <summary>
        /// Imputes missing intensity value for each protein
        /// </summary>
        private void ImputeData(ProteinRowInfo proteinRowInfo, double[] samplesMeanIntensityValue, double[] samplesStandardDeviation,
            List<string> samplesFileNames, double[] missingFactor, int[] numberOfIntensityValuesInSample, double meanFraction)
        {
            Dictionary<string, double> samplesintensityData = proteinRowInfo.getSamplesIntensityValues();

            for (int i = 0; i < samplesFileNames.Count; i++)
            {
                if (samplesintensityData[samplesFileNames[i]] == 0)
                {

                    double imputedFraction = missingFactor[i] / (numberOfIntensityValuesInSample[i] + missingFactor[i]);
                    if (imputedFraction > 0.5) continue;
                    double imputedProbability = imputedFraction / (1 - imputedFraction);
                    double standardDeviationFraction = Math.Max(2 * imputedFraction, 0.3);
                    double stdDevFraction = 0.6 * (1 - (imputedFraction * imputedFraction));
                    MathNet.Numerics.Distributions.Normal probabilityDist = new Normal(samplesMeanIntensityValue[i], standardDeviationFraction);
                    double probabilitySetPoint = probabilityDist.Density(samplesMeanIntensityValue[i] + stdDevFraction * standardDeviationFraction);
                    double yCoordinate = imputedProbability * probabilitySetPoint;
                    double deltaX = standardDeviationFraction * stdDevFraction;
                    MathNet.Numerics.Distributions.Normal xCoord = new Normal(samplesMeanIntensityValue[i], samplesStandardDeviation[i]);
                    double deltaMu = xCoord.InverseCumulativeDistribution(yCoordinate);
                    double meanDownshift = (deltaMu - deltaX * meanFraction);


                    MathNet.Numerics.Distributions.Normal normalDist = new Normal(meanDownshift, standardDeviationFraction);
                    double imputeVal = normalDist.Sample();
                    samplesintensityData[samplesFileNames[i]] = imputeVal;
                }
            }
        }

        /// <summary>
        /// Computes all required data and then imputes missing intensity values of protein
        /// </summary>
        public void RunImputationProcess(List<ProteinRowInfo> allProteinInfo, List<string> samplesFileNames, double meanFraction)
        {

            double[] samplesAllIntensityValuesSum = new double[samplesFileNames.Count];
            double[] samplesMeanIntensityValue = new double[samplesFileNames.Count];
            double[] samplesStandardDeviationNumerators = new double[samplesFileNames.Count];
            double[] samplesStandardDeviation = new double[samplesFileNames.Count];
            double[] missingFactor = new double[samplesFileNames.Count]; // MIGHT HAVE TO REMOVE
            int[] numberOfIntensityValuesInSample = new int[samplesFileNames.Count];

            for (int i = 0; i < allProteinInfo.Count; i++)
            {
                ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                CalculateNumberAndSumOfIntensityValuesOfSample(proteinRowInfo, samplesAllIntensityValuesSum, samplesFileNames,
                    numberOfIntensityValuesInSample);
            }

            for (int i = 0; i < samplesAllIntensityValuesSum.Length; i++)
            {
                // computes the mean of intensity values for each sample 
                samplesMeanIntensityValue[i] = samplesAllIntensityValuesSum[i] / numberOfIntensityValuesInSample[i];
            }

            for (int i = 0; i < allProteinInfo.Count; i++)
            {
                ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                CalculateSampleStandardDeviationNumerator(proteinRowInfo, samplesStandardDeviationNumerators, samplesFileNames,
                    samplesMeanIntensityValue);
            }

            for (int i = 0; i < allProteinInfo.Count; i++)
            {
                ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                CalculateMissingIntensityValuesInSample(proteinRowInfo, missingFactor, samplesFileNames);
            }

            for (int i = 0; i < samplesStandardDeviationNumerators.Length; i++)
            {
                // computes the standard deviation of intensity values for each sample 
                samplesStandardDeviation[i] = Math.Sqrt(samplesStandardDeviationNumerators[i] / (allProteinInfo.Count - missingFactor[i]));
            }

            for (int i = 0; i < allProteinInfo.Count; i++)
            {
                ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                ImputeData(proteinRowInfo, samplesMeanIntensityValue, samplesStandardDeviation, samplesFileNames, missingFactor,
                    numberOfIntensityValuesInSample, meanFraction);
            }
        }
    }
}
