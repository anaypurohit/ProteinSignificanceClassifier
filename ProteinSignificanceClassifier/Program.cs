using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ProteinSignificanceClassifier
{
    class Program
    {
        /// <summary>
        /// Parses the Experimental Design files outputted from FlashLFQ
        /// Determines the conditions and their corresponding samples
        /// </summary>
        private Dictionary<string, List<string>> ExpermientalDesignParser(string experimentalDesignLocation)
        {
            Dictionary<string, List<string>> samplefileConditionRelation = new Dictionary<string, List<string>>();
            int skipFirstLine = 0;
            foreach (var line in File.ReadLines(experimentalDesignLocation))
            {
                if (skipFirstLine == 0)
                {
                    skipFirstLine = 1;
                    continue;
                }
                var experimentInfo = line.Split('\t');
                string sampleFileName = experimentInfo[0];
                string sampleCondition = experimentInfo[1];

                samplefileConditionRelation.TryGetValue(sampleCondition, out List<string> sampleFiles);
                if (sampleFiles == null)
                {
                    samplefileConditionRelation.Add(sampleCondition, new List<string>() { sampleFileName });
                }
                else
                {
                    sampleFiles.Add(sampleFileName);
                }

            }
            return samplefileConditionRelation;
        }


        /// <summary>
        /// Parses the QuantifiedProteins files outputted from FlashLFQ
        /// Filters out proteins with >="numNans" number of missing Intensity Values 
        /// </summary>
        private void ProteinDataParser(List<ProteinRowInfo> allProteinInfo, int maxInvalidIntensityValues,
            List<string> samplesFileNames, string quantifiedProteinFileLocation)
        {
            int startIntensitycolNum = 0;
            int endIntensitycolNum = 0;
            int readFirstRow = 0;
            int maxAllowedInvalidIntensityValues = maxInvalidIntensityValues;

            foreach (var line in File.ReadLines(quantifiedProteinFileLocation))
            {
                var proteinRow = line.Split('\t');

                // determines location of Intensity values and stores samplefileNames based on first line in QuantifiedProteins file
                if (readFirstRow == 0)
                {
                    for (int i = 0; i < proteinRow.Length; i++)
                    {
                        if (proteinRow[i].StartsWith("Intensity"))
                        {
                            samplesFileNames.Add(proteinRow[i]);
                            if (startIntensitycolNum == 0) startIntensitycolNum = i;
                            endIntensitycolNum = i;
                        }
                    }
                    readFirstRow = 1;
                    continue;
                }

                // determine how many invalid Intensity values the given protein has
                int invalidIntensityValuesForProtein = 0;
                for (int i = startIntensitycolNum; i <= endIntensitycolNum; i++)
                {
                    if (string.IsNullOrWhiteSpace(proteinRow[i]) || Convert.ToDouble(proteinRow[i]) == 0) invalidIntensityValuesForProtein++;
                }
                if (invalidIntensityValuesForProtein >= maxAllowedInvalidIntensityValues) continue;

                // we know now that sample is valid so will be stored
                ProteinRowInfo proteinRowInfo = new ProteinRowInfo();
                proteinRowInfo.setProteinID(proteinRow[1]); // set protein ID

                int sampleFileNumberTracker = 0; // used to add the correct intensity value
                for (int i = startIntensitycolNum; i <= endIntensitycolNum; i++)
                {
                    if (string.IsNullOrWhiteSpace(proteinRow[i]))
                    {
                        proteinRowInfo.setSamplesIntensityValues(samplesFileNames[sampleFileNumberTracker], 0);
                    }
                    else
                    {
                        if (Convert.ToDouble(proteinRow[i]) == 0)
                        {
                            proteinRowInfo.setSamplesIntensityValues(samplesFileNames[sampleFileNumberTracker], 0);
                        }
                        else
                        {
                            proteinRowInfo.setSamplesIntensityValues(samplesFileNames[sampleFileNumberTracker], Math.Log(Convert.ToDouble(proteinRow[i]), 2));
                        }
                    }
                    sampleFileNumberTracker++;
                }

                //proper row sample added to list
                allProteinInfo.Add(proteinRowInfo);
            }
        }


        /// <summary>
        /// Used to resize the permuted N values containing array so that its size is equal 
        /// to the original N values array size
        /// </summary>
        private void resizePermutedArray(List<double> permutedNValues, int sizeDifference)
        {
            if (sizeDifference <= 0)
            {
                return;
            }
            else
            {
                permutedNValues.Sort();
                var trackElements = (Convert.ToDouble(permutedNValues.Count) - 1) / sizeDifference;
                double loopMaintainer = trackElements;


                while (trackElements < permutedNValues.Count)
                {
                    double roundedNumber = Math.Round(trackElements);
                    int roundedInt = Convert.ToInt32(roundedNumber);
                    permutedNValues[roundedInt] = -1;
                    trackElements = trackElements + loopMaintainer;
                }

                permutedNValues.RemoveAll(value => value == -1);
            }
        }

        /// <summary>
        /// This method is used to generate all possible combinations of choosing 2 conditions from all possible conditions
        /// </summary>
        private List<List<string>> GenerateAllCombinationsOfTwoConditions(List<string> conditions)
        {
            List<List<string>> allTwoConditionCombinations = new List<List<string>>();
            for (int i = 0; i < conditions.Count; i++)
            {
                for (int j = 0; j < conditions.Count; j++)
                {
                    if (i < j)
                    {
                        allTwoConditionCombinations.Add(new List<string> { conditions[i], conditions[j] });
                    }
                }
            }
            return allTwoConditionCombinations;
        }


        /// <summary>
        /// Determines which proteins are significant based on N Value Threshold and prints out the classifications
        /// </summary>
        private void printSignificantProtein(List<ProteinRowInfo> allProteinInfo, List<double> actualNValues, double nValueThreshold,
            string printSignificantProteinsLocation)
        {

            //"C:/Users/Anay/Desktop/UW Madison/Smith Lab/Project 1/ConsoleApp1/ChangedSemicolon11.csv"
            using (StreamWriter writetext = new StreamWriter(printSignificantProteinsLocation))
            {
                for (int i = 0; i < actualNValues.Count; i++)
                {
                    ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                    if (actualNValues[i] < nValueThreshold)
                    {
                        writetext.Write(proteinRowInfo.getProteinID() + ", " + "Not Significant");
                        writetext.WriteLine();
                    }
                    else
                    {
                        writetext.Write(proteinRowInfo.getProteinID() + ", " + "Significant");
                        writetext.WriteLine();
                    }
                }
            }
        }

        public static void Main(string[] args)
        {
            double sOValue = 0.1;//0.3;
            double meanFraction = 0.1;//5.7; //2.4
            int maxSignificantCount = 0;

            while (meanFraction < 6)
            {
                while (sOValue < 2)
                {
                    for (int k = 0; k < 8; k++)
                    {
                        Program dataManipulate = new Program();

                        // Parse the ExperimentalDesign File to get info of samples and conditions they belong to
                        Dictionary<string, List<string>> samplefileConditionRelation = dataManipulate.ExpermientalDesignParser("C:/Users/Anay/Desktop/UW Madison/Smith Lab/Spectra Data/ExperimentalDesign.tsv");

                        List<ProteinRowInfo> allProteinInfo = new List<ProteinRowInfo>();
                        List<string> samplesFileNames = new List<string>();
                        int maxInvalidIntensityValues = k;
                        dataManipulate.ProteinDataParser(allProteinInfo, maxInvalidIntensityValues,
                            samplesFileNames, "C:/Users/Anay/Desktop/UW Madison/Smith Lab/Spectra Data/FlashLFQ_2020-04-26-17-39-35/QuantifiedProteins.tsv");

                        // imputes missing intensity values for each protein
                        ImputationProcess imputationProcess = new ImputationProcess();
                        imputationProcess.RunImputationProcess(allProteinInfo, samplesFileNames, meanFraction);

                        List<double> actualNValues = new List<double>(); // will store actual(real) N values
                        List<double> permutedNValues = new List<double>(); // will store permuted(fake) N values
                        StatisticalTests statisticalTests = new StatisticalTests();

                        // Compute actual and permuted N Values for each protein by considering every condition pair
                        for (int i = 0; i < allProteinInfo.Count; i++)
                        {
                            ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                            Dictionary<string, double> samplesintensityData = proteinRowInfo.getSamplesIntensityValues();

                            // get all conditions and pair them up for the T Tests and Permutation Testing
                            List<string> allConditions = new List<string>(samplefileConditionRelation.Keys);
                            List<List<string>> allTwoConditionCombinations = dataManipulate.GenerateAllCombinationsOfTwoConditions(allConditions);

                            foreach (List<string> conditionPair in allTwoConditionCombinations)
                            {
                                string firstCondition = conditionPair[0];
                                List<string> firstConditionAssociatedSamples = samplefileConditionRelation.GetValueOrDefault(firstCondition);
                                string secondCondition = conditionPair[1];
                                List<string> secondConditionAssociatedSamples = samplefileConditionRelation.GetValueOrDefault(secondCondition);
                                List<double> proteinFirstConditionIntensityValues = new List<double>();
                                List<double> proteinSecondConditionIntensityValues = new List<double>();

                                // get the protein's intensity values corresponding to the chosen pair of conditions
                                foreach (string sampleFileName in samplesFileNames)
                                {
                                    if (firstConditionAssociatedSamples.Contains(sampleFileName))
                                    {
                                        proteinFirstConditionIntensityValues.Add(samplesintensityData[sampleFileName]);
                                    }
                                    if (secondConditionAssociatedSamples.Contains(sampleFileName))
                                    {
                                        proteinSecondConditionIntensityValues.Add(samplesintensityData[sampleFileName]);
                                    }
                                }

                                // Compute actual(real) N Values with the chosen pair of conditions using T-Tests and
                                // store in actualNValues array
                                statisticalTests.GetNValueUsingTTest(proteinFirstConditionIntensityValues, proteinSecondConditionIntensityValues,
                                    actualNValues, sOValue);

                                // Compute permuted(fake) N Values with the chosen pair of conditions using T-Tests and 
                                // store in permutedNValues array
                                statisticalTests.GetNValueUsingPermutationtests(proteinFirstConditionIntensityValues, proteinSecondConditionIntensityValues,
                                    permutedNValues, sOValue);
                            }
                        }

                        // makes the permuted N values list and the actual N Values list of the same size
                        dataManipulate.resizePermutedArray(permutedNValues, permutedNValues.Count() - actualNValues.Count());

                        // Copy of the actual N values which will be used when determind the N Value threshold for target FDR 
                        List<double> actualNValuesCopy = new List<double>();
                        for (int i = 0; i < actualNValues.Count; i++)
                        {
                            actualNValuesCopy.Add(actualNValues[i]);
                        }
                        // get the threshold at which we will filter out the significant proteins
                        double nValueThreshold = statisticalTests.calculateNvaluethreshold(actualNValuesCopy, permutedNValues, 0.05);

                        // determine number of signifcant proteins detected
                        int newSignificantCount = actualNValues.Count(x => x >= nValueThreshold);
                        if (newSignificantCount > maxSignificantCount)
                        {
                            maxSignificantCount = newSignificantCount;
                        }
                    }
                    sOValue = sOValue + 0.1;
                }
                sOValue = 0.1;
                meanFraction = meanFraction + 0.3;
            }
        }
    }
}
