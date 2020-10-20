using Accord.Statistics.Testing;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ProteinSignificanceClassifier
{
    /// <summary>
    /// This class contains the statsitical tests used in the program
    /// </summary>
    public class StatisticalTests
    {
        /// <summary>
        /// This method is used to generate all possible combinations of choosing 2 indicies from
        /// the given indicies
        /// </summary>
        private static List<List<int>> GenerateAllCombinationsOfTwoIndices(int[] indices)
        {
            List<List<int>> allTwoIndicesCombination = new List<List<int>>();
            for (int i = 0; i < indices.Length; i++)
            {
                for (int j = 0; j < indices.Length; j++)
                {
                    if (i < j)
                    {
                        allTwoIndicesCombination.Add(new List<int> { i, j });
                    }
                }
            }
            return allTwoIndicesCombination;
        }

        /// <summary>
        /// Calculates the standard deviation of intensity values of the protein
        /// </summary>
        private double CalculateProteinIntensityValuesStandardDeviation(double[] vals, double mean)
        {
            double stdDev = 0;
            for (int i = 0; i < vals.Length; i++)
            {
                stdDev = stdDev + (vals[i] - mean) * (vals[i] - mean);
            }
            stdDev = stdDev / vals.Length;
            stdDev = Math.Sqrt(stdDev);
            return stdDev;

        }

        /// <summary>
        /// Calculates the mean of intensity values of the protein
        /// </summary>
        private double CalculateProteinMeanIntensityValue(double[] vals)
        {
            double mean = 0;
            for (int i = 0; i < vals.Length; i++)
            {
                mean = mean + vals[i];
            }
            mean = mean / vals.Length;
            return mean;
        }


        /// <summary>
        /// Calculates the N Value for each protein based on its intensity values in two conditions
        /// </summary>
        public void GetNValueUsingTTest(List<double> proteinFirstConditionIntensityValues, List<double> proteinSecondConditionIntensityValues,
            List<double> actualNValues, List<double> actualPValues, List<double> actualLogFoldChange, double sOValue)
        {
            double[] firstConditionIntensityValues = new double[proteinFirstConditionIntensityValues.Count];
            double[] secondConditionIntensityValues = new double[proteinSecondConditionIntensityValues.Count];

            for (int i = 0; i < proteinFirstConditionIntensityValues.Count; i++)
            {
                firstConditionIntensityValues[i] = proteinFirstConditionIntensityValues[i];
            }

            for (int i = 0; i < proteinSecondConditionIntensityValues.Count; i++)
            {
                secondConditionIntensityValues[i] = proteinSecondConditionIntensityValues[i];
            }


            double firstConditionIntensityMean = CalculateProteinMeanIntensityValue(firstConditionIntensityValues);
            double secondConditionIntensityMean = CalculateProteinMeanIntensityValue(secondConditionIntensityValues);
            double firstConditionIntensityStandardDev = CalculateProteinIntensityValuesStandardDeviation(firstConditionIntensityValues, firstConditionIntensityMean);
            double secondConditionIntensityStandardDev = CalculateProteinIntensityValuesStandardDeviation(secondConditionIntensityValues, secondConditionIntensityMean);
            double firstConditionIntensityVariance = firstConditionIntensityStandardDev * firstConditionIntensityStandardDev;
            double secondConditionIntensityVariance = secondConditionIntensityStandardDev * secondConditionIntensityStandardDev;

            var ftest = new FTest(firstConditionIntensityVariance, secondConditionIntensityVariance, proteinFirstConditionIntensityValues.Count - 1, proteinSecondConditionIntensityValues.Count - 1);
            bool significant = ftest.Significant; // gets whether null hypothesis can be rejected

            // Create two tailed t test to get p values
            TwoSampleTTest ttest = new TwoSampleTTest(firstConditionIntensityValues, secondConditionIntensityValues, !significant);
            double pValue = ttest.PValue;
            double logpValue = -Math.Log10(pValue);
            double logfoldChange = secondConditionIntensityMean - firstConditionIntensityMean;
            double nValue = (logpValue * (logfoldChange * logfoldChange - sOValue * sOValue)) / ((logfoldChange) * (logfoldChange));

            actualNValues.Add(nValue);
            actualLogFoldChange.Add(logfoldChange);
            actualPValues.Add(pValue);
        }

        /// <summary>
        /// Generates fake N Value's for each protein based on its intensity values in two conditions using Permutation Testing
        /// </summary>
        public void GetNValueUsingPermutationtests(List<double> proteinFirstConditionIntensityValues, List<double> proteinSecondConditionIntensityValues,
             List<double> permutedNValues, double sOValue)
        {
            double[] firstConditionIntensityValues = new double[proteinFirstConditionIntensityValues.Count];
            double[] secondConditionIntensityValues = new double[proteinSecondConditionIntensityValues.Count];

            for (int i = 0; i < proteinFirstConditionIntensityValues.Count; i++)
            {
                firstConditionIntensityValues[i] = proteinFirstConditionIntensityValues[i];
            }

            for (int i = 0; i < proteinSecondConditionIntensityValues.Count; i++)
            {
                secondConditionIntensityValues[i] = proteinSecondConditionIntensityValues[i];
            }


            int[] indicesOfFirstConditionIntensityValues = new int[firstConditionIntensityValues.Length];
            int[] indicesOfSecondConditionIntensityValues = new int[secondConditionIntensityValues.Length];
            for (int i = 0; i < proteinFirstConditionIntensityValues.Count; i++)
            {
                indicesOfFirstConditionIntensityValues[i] = i;
            }
            for (int i = 0; i < proteinSecondConditionIntensityValues.Count; i++)
            {
                indicesOfSecondConditionIntensityValues[i] = i;
            }

            List<List<int>> allTwoIndiciesCombinationsFromFirstCondition = GenerateAllCombinationsOfTwoIndices(indicesOfFirstConditionIntensityValues);
            List<List<int>> allTwoIndiciesCombinationsFromSecondCondition = GenerateAllCombinationsOfTwoIndices(indicesOfSecondConditionIntensityValues);

            int count = 0;
            foreach (var twoIndiciesCombinationEntryFromFirstCondition in allTwoIndiciesCombinationsFromFirstCondition)
            {
                foreach (var twoIndiciesCombinationEntryFromSecondCondition in allTwoIndiciesCombinationsFromSecondCondition)
                {
                    // these the new arrays which will be made after swapping intensity values between the two conditions
                    double[] swappedFirstConditionIntensityValues = new double[firstConditionIntensityValues.Length];
                    double[] swappedSecondConditionIntensityValues = new double[secondConditionIntensityValues.Length];
                    int swappedFirstConditionArrayTracker = 0;
                    int swappedSecondConditionArrayTracker = 0;

                    int[] indiciesToSwapFromFirstCondition = new int[2];
                    int[] indiciesToSwapFromSecondCondition = new int[2];
                    int removeIndiciesFirstConditionTracker = 0;
                    int removeIndiciesSecondConditionTracker = 0;

                    // store the indices, corresponding to intensity values, to be swapped from first condition
                    foreach (var index in twoIndiciesCombinationEntryFromFirstCondition)
                    {
                        indiciesToSwapFromFirstCondition[removeIndiciesFirstConditionTracker] = index;
                        removeIndiciesFirstConditionTracker++;
                    }

                    // store the indices, corresponding to intensity values, to be swapped from second condition
                    foreach (var index in twoIndiciesCombinationEntryFromSecondCondition)
                    {
                        indiciesToSwapFromSecondCondition[removeIndiciesSecondConditionTracker] = index;
                        removeIndiciesSecondConditionTracker++;
                    }

                    // add the intensity values to be swapped from first condition the second condition
                    for (int j = 0; j < indiciesToSwapFromFirstCondition.Count(); j++)
                    {
                        swappedSecondConditionIntensityValues[swappedSecondConditionArrayTracker] = firstConditionIntensityValues[indiciesToSwapFromFirstCondition[j]];
                        swappedSecondConditionArrayTracker++;
                    }


                    // add the intensity values to be swapped from second condition the first condition
                    for (int j = 0; j < indiciesToSwapFromSecondCondition.Count(); j++)
                    {
                        swappedFirstConditionIntensityValues[swappedFirstConditionArrayTracker] = secondConditionIntensityValues[indiciesToSwapFromSecondCondition[j]];
                        swappedFirstConditionArrayTracker++;
                    }

                    // now we add the remaining intensity values from the first condition to the swappedFirstCondition Array
                    for (int j = 0; j < firstConditionIntensityValues.Count(); j++)
                    {
                        if (indiciesToSwapFromFirstCondition.Contains(j)) continue;
                        swappedFirstConditionIntensityValues[swappedFirstConditionArrayTracker] = firstConditionIntensityValues[j];
                        swappedFirstConditionArrayTracker++;
                    }

                    // now we add the remaining intensity values from the second condition to the swappedSecondCondition Array
                    for (int j = 0; j < secondConditionIntensityValues.Count(); j++)
                    {
                        if (indiciesToSwapFromSecondCondition.Contains(j)) continue;
                        swappedSecondConditionIntensityValues[swappedSecondConditionArrayTracker] = secondConditionIntensityValues[j];
                        swappedSecondConditionArrayTracker++;
                    }


                    // at this stage we have the newly made swapped arrays with mixture of groups.
                    // need to proceed with T tests for these groups to generate fake p values.

                    double firstConditionIntensityMean = CalculateProteinMeanIntensityValue(swappedFirstConditionIntensityValues);
                    double secondConditionIntensityMean = CalculateProteinMeanIntensityValue(swappedSecondConditionIntensityValues);
                    double firstConditionIntensityStandardDev = CalculateProteinIntensityValuesStandardDeviation(swappedFirstConditionIntensityValues, firstConditionIntensityMean);
                    double secondConditionIntensityStandardDev = CalculateProteinIntensityValuesStandardDeviation(swappedSecondConditionIntensityValues, secondConditionIntensityMean);
                    double firstConditionIntensityVariance = firstConditionIntensityStandardDev * firstConditionIntensityStandardDev;
                    double secondConditionIntensityVariance = secondConditionIntensityStandardDev * secondConditionIntensityStandardDev;

                    var ftest = new FTest(firstConditionIntensityVariance, secondConditionIntensityVariance, swappedFirstConditionIntensityValues.Length - 1, swappedSecondConditionIntensityValues.Length - 1);
                    bool significant = ftest.Significant; // gets whether null hypothesis can be rejected

                    // Create two tailed t test to get p values
                    TwoSampleTTest ttest = new TwoSampleTTest(swappedFirstConditionIntensityValues, swappedSecondConditionIntensityValues, !significant);
                    double pValue = ttest.PValue;
                    double logpValue = -Math.Log10(pValue);
                    double logfoldChange = secondConditionIntensityMean - firstConditionIntensityMean;

                    permutedNValues.Add((logpValue * (logfoldChange * logfoldChange - sOValue * sOValue)) / ((logfoldChange) * (logfoldChange)));
                }
                count++;
                if (count == 2) break;
            }
        }

        /// <summary>
        /// Used to Determine the threshold N value over which a protein would be classified as being significant.
        /// It is used to determine the N value threhold for FDR purposes. It does so combining the permuted N values 
        /// and actual N values in descending order and determines when (perumted N values)/(actual N Values) > Target FDR.
        /// The N Value when that happens becomes the threshold for achieving the target FDR.
        /// </summary>
        public double calculateNvaluethreshold(List<double> actualNValues, List<double> perumtedNValues, double FDR)
        {
            actualNValues.Sort();
            actualNValues.Reverse();
            perumtedNValues.Sort();
            perumtedNValues.Reverse();
            double countActualNVals = 0;
            double countPermutedNVals = 0;
            int actualNvalsArrayTracker = 0;
            int permutedNvalsArrayTracker = 0;
            List<double> combinedNVals = new List<double>();
            while (true)
            {
                if (actualNvalsArrayTracker <= actualNValues.Count() - 1 && permutedNvalsArrayTracker <= perumtedNValues.Count() - 1)
                {
                    if (actualNValues[actualNvalsArrayTracker] > perumtedNValues[permutedNvalsArrayTracker])
                    {
                        combinedNVals.Add(actualNValues[actualNvalsArrayTracker]);
                        countActualNVals++;
                        actualNvalsArrayTracker++;
                    }
                    else if (actualNValues[actualNvalsArrayTracker] < perumtedNValues[permutedNvalsArrayTracker])
                    {
                        combinedNVals.Add(perumtedNValues[permutedNvalsArrayTracker]);
                        countPermutedNVals++;
                        permutedNvalsArrayTracker++;
                        if (countActualNVals == 0 || ((countPermutedNVals / countActualNVals) > FDR))
                        {
                            return combinedNVals[combinedNVals.Count - 1];
                        }
                    }
                    else
                    {
                        combinedNVals.Add(actualNValues[actualNvalsArrayTracker]);
                        combinedNVals.Add(perumtedNValues[permutedNvalsArrayTracker]);
                        countActualNVals++;
                        countPermutedNVals++;
                        actualNvalsArrayTracker++;
                        permutedNvalsArrayTracker++;
                        if ((countPermutedNVals / countActualNVals) > FDR)
                        {
                            return combinedNVals[combinedNVals.Count - 1];
                        }
                    }
                }
                else if (actualNvalsArrayTracker > actualNValues.Count() - 1 && permutedNvalsArrayTracker <= perumtedNValues.Count() - 1)
                {
                    combinedNVals.Add(perumtedNValues[permutedNvalsArrayTracker]);
                    countPermutedNVals++;
                    permutedNvalsArrayTracker++;
                    if (countActualNVals == 0 || (countPermutedNVals / countActualNVals) > FDR)
                    {
                        return combinedNVals[combinedNVals.Count - 1];
                    }
                }
                else if (actualNvalsArrayTracker <= actualNValues.Count() - 1 && permutedNvalsArrayTracker > perumtedNValues.Count() - 1)
                {
                    combinedNVals.Add(actualNValues[actualNvalsArrayTracker]);
                    countActualNVals++;
                    actualNvalsArrayTracker++;
                }
                else
                {
                    // need to return last element in combined array as ratio never exceeded
                    return combinedNVals[combinedNVals.Count() - 1];
                }
            }
        }
    }
}
