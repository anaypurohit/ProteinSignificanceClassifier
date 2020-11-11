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
                for (int j = i + 1; j < indices.Length; j++)
                {
                    allTwoIndicesCombination.Add(new List<int> { i, j });
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
            for (int i = 0; i < vals.Length; i++) //better to use foreach here, because you have to lookup vals[i] twice.
            {
                stdDev += (vals[i] - mean) * (vals[i] - mean);
            }
            stdDev /= vals.Length;
            stdDev = Math.Sqrt(stdDev);

            return stdDev;
        }

        /// <summary>
        /// Calculates the mean of intensity values of the protein
        /// </summary>
        private double CalculateProteinMeanIntensityValue(double[] vals)
        {
            //Delete this method. Just use vals.Average();
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
        /// //need explanation on what an N value is
        public void GetNValueUsingTTest(List<double> proteinFirstConditionIntensityValues, List<double> proteinSecondConditionIntensityValues,
            List<double> actualNValues, List<double> actualPValues, List<double> actualLogFoldChange, double sOValue)
        {
            //why are these being cloned?
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

            //sometimes less is more with variable names. For example, "mean1" vs "firstConditionIntensityMean" improves readability without diminishing clarity
            double firstConditionIntensityMean = CalculateProteinMeanIntensityValue(firstConditionIntensityValues);
            double secondConditionIntensityMean = CalculateProteinMeanIntensityValue(secondConditionIntensityValues);
            double firstConditionIntensityStandardDev = CalculateProteinIntensityValuesStandardDeviation(firstConditionIntensityValues, firstConditionIntensityMean);
            double secondConditionIntensityStandardDev = CalculateProteinIntensityValuesStandardDeviation(secondConditionIntensityValues, secondConditionIntensityMean);
            double firstConditionIntensityVariance = firstConditionIntensityStandardDev * firstConditionIntensityStandardDev;
            double secondConditionIntensityVariance = secondConditionIntensityStandardDev * secondConditionIntensityStandardDev;

            //don't need to save "ftest" if you're never going to use it again
            //write a note here what the purpose of the f-test is. It's not as well known as the t-test
            var ftest = new FTest(firstConditionIntensityVariance, secondConditionIntensityVariance, proteinFirstConditionIntensityValues.Count - 1, proteinSecondConditionIntensityValues.Count - 1);
            bool significant = ftest.Significant; // gets whether null hypothesis can be rejected

            // Create two tailed t test to get p values
            TwoSampleTTest ttest = new TwoSampleTTest(firstConditionIntensityValues, secondConditionIntensityValues, !significant);
            double pValue = ttest.PValue;
            double logpValue = -Math.Log10(pValue);
            double logfoldChange = secondConditionIntensityMean - firstConditionIntensityMean;
            //need a note on where this came from
            double nValue = (logpValue * (logfoldChange * logfoldChange - sOValue * sOValue)) / ((logfoldChange) * (logfoldChange));

            //these are interesting names that imply the existence of fake nvalues, foldchanges, and pvalues.
            //how about "allNValues"?
            actualNValues.Add(nValue);
            actualLogFoldChange.Add(logfoldChange);
            actualPValues.Add(pValue);
        }

        /// <summary>
        /// Generates fake N Value's for each protein based on its intensity values in two conditions using Permutation Testing
        /// </summary>
        /// //alright, you win, I guess there are fake nvalues. Maybe a better nomenclature would be observedNValues vs permutedNValues?
        /// There's a lot of code duplication between this method and GetNValueUsingTTest. Try creating permutations and then calling GetNValueUsingTTest for each one.
        public void GetNValueUsingPermutationtests(List<double> proteinFirstConditionIntensityValues, List<double> proteinSecondConditionIntensityValues,
             List<double> permutedNValues, double sOValue)
        {
            //why are these being cloned?
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
            actualNValues = actualNValues.OrderByDescending(x => x).ToList();
            perumtedNValues = perumtedNValues.OrderByDescending(x => x).ToList();
            double countActualNVals = 0;
            double countPermutedNVals = 0; //much nicer than fake :)
            int actualNValsArrayTracker = 0;
            int permutedNValsArrayTracker = 0;
            List<double> combinedNVals = new List<double>();
            while (true) //gross. Can we refactor to remove this?
            {
                if (actualNValsArrayTracker <= actualNValues.Count - 1 && permutedNValsArrayTracker <= perumtedNValues.Count - 1)
                {
                    if (actualNValues[actualNValsArrayTracker] > perumtedNValues[permutedNValsArrayTracker])
                    {
                        combinedNVals.Add(actualNValues[actualNValsArrayTracker]);
                        countActualNVals++;
                        actualNValsArrayTracker++;
                    }
                    else if (actualNValues[actualNValsArrayTracker] < perumtedNValues[permutedNValsArrayTracker])
                    {
                        combinedNVals.Add(perumtedNValues[permutedNValsArrayTracker]);
                        countPermutedNVals++;
                        permutedNValsArrayTracker++;
                        if (countActualNVals == 0 || ((countPermutedNVals / countActualNVals) > FDR))
                        {
                            return combinedNVals[combinedNVals.Count - 1];
                        }
                    }
                    else
                    {
                        combinedNVals.Add(actualNValues[actualNValsArrayTracker]);
                        combinedNVals.Add(perumtedNValues[permutedNValsArrayTracker]);
                        countActualNVals++;
                        countPermutedNVals++;
                        actualNValsArrayTracker++;
                        permutedNValsArrayTracker++;
                        if ((countPermutedNVals / countActualNVals) > FDR)
                        {
                            return combinedNVals[combinedNVals.Count - 1];
                        }
                    }
                }
                else if (actualNValsArrayTracker > actualNValues.Count - 1 && permutedNValsArrayTracker <= perumtedNValues.Count - 1)
                {
                    combinedNVals.Add(perumtedNValues[permutedNValsArrayTracker]);
                    countPermutedNVals++;
                    permutedNValsArrayTracker++;
                    if (countActualNVals == 0 || (countPermutedNVals / countActualNVals) > FDR)
                    {
                        return combinedNVals[combinedNVals.Count - 1];
                    }
                }
                else if (actualNValsArrayTracker <= actualNValues.Count - 1 && permutedNValsArrayTracker > perumtedNValues.Count - 1)
                {
                    combinedNVals.Add(actualNValues[actualNValsArrayTracker]);
                    countActualNVals++;
                    actualNValsArrayTracker++;
                }
                else
                {
                    // need to return last element in combined array as ratio never exceeded
                    return combinedNVals[combinedNVals.Count - 1];
                }
            }
        }
    }
}
