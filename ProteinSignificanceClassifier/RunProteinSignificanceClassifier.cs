using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace ProteinSignificanceClassifier
{
    public class RunProteinSignificanceClassifier
    {
        /// <summary>
        /// Parses the Experimental Design files outputted from FlashLFQ
        /// Determines the conditions and their corresponding samples
        /// </summary>
        public Dictionary<string, List<string>> ExpermientalDesignParser(string experimentalDesignLocation, ref int numberSamples)
        {
            Dictionary<string, List<string>> samplefileConditionRelation = new Dictionary<string, List<string>>();

            //read in experimental design
            string[] lines = File.ReadLines(experimentalDesignLocation).ToArray();

            for (int l = 1; l < lines.Length; l++)
            {
                var experimentInfo = lines[l].Split('\t');
                string sampleFileName = experimentInfo[0];
                string sampleCondition = experimentInfo[1];

                samplefileConditionRelation.TryGetValue(sampleCondition, out List<string> sampleFiles);
                if (sampleFiles == null)
                {
                    samplefileConditionRelation.Add(sampleCondition, new List<string>() { "Intensity_" + sampleFileName });
                }
                else
                {
                    sampleFiles.Add("Intensity_" + sampleFileName);
                }
                numberSamples++;
            }
            return samplefileConditionRelation;
        }


        /// <summary>
        /// Parses the QuantifiedProteins files outputted from FlashLFQ
        /// Filters out proteins with >="numNans" number of missing Intensity Values and generates a dictionary file
        /// mapping all conditions to their sampleFileName
        /// </summary>
        public void ProteinDataParser(List<ProteinRowInfo> allProteinInfo, int maxAllowedMissingValues,
            List<string> samplesFileNames, string quantifiedProteinFileLocation)
        {
            int startIntensityIndex = 0;
            int endIntensityIndex = 0;
            int proteinIDIndex = 0;

            //read in results
            string[] lines = File.ReadLines(quantifiedProteinFileLocation).ToArray();

            //parse header
            string[] header = lines[0].Split('\t');
            for (int i = 0; i < header.Length; i++)
            {
                if (header[i].StartsWith("Intensity"))
                {
                    samplesFileNames.Add(header[i]);
                    if (startIntensityIndex == 0)
                    {
                        startIntensityIndex = i;
                    }
                    endIntensityIndex = i;
                }

                if (header[i].StartsWith("Protein"))
                {
                    proteinIDIndex = i;
                }
            }

            //parse data
            for (int l = 1; l < lines.Length; l++)
            {
                var proteinRow = lines[l].Split('\t');

                // determine how many missing Intensity values the given protein has
                int numMissingValues = 0;
                for (int i = startIntensityIndex; i <= endIntensityIndex; i++)
                {
                    if (string.IsNullOrWhiteSpace(proteinRow[i]) || Convert.ToDouble(proteinRow[i]) == 0)
                    {
                        numMissingValues++;
                    }
                }

                //if there are enough valid values
                if (numMissingValues <= maxAllowedMissingValues)
                {
                    // we know now that sample is valid so will be stored
                    ProteinRowInfo proteinRowInfo = new ProteinRowInfo();
                    proteinRowInfo.ProteinID = proteinRow[proteinIDIndex]; // set protein ID

                    int sampleFileNumberTracker = 0; // used to add the correct intensity value
                    for (int i = startIntensityIndex; i <= endIntensityIndex; i++)
                    {
                        string writtenValue = proteinRow[i];
                        double intensity = string.IsNullOrWhiteSpace(writtenValue) || Convert.ToDouble(writtenValue) == 0 ? 
                            0 : Math.Log(Convert.ToDouble(writtenValue), 2);
                        proteinRowInfo.SamplesIntensityData[samplesFileNames[sampleFileNumberTracker]] = intensity;
                       
                        sampleFileNumberTracker++;
                    }

                    //proper row sample added to list
                    allProteinInfo.Add(proteinRowInfo);
                }
            }
        }


        /// <summary>
        /// Used to resize the permuted N values containing array so that its size is equal 
        /// to the original N values array size
        /// </summary>
        public void ResizePermutedArray(List<double> permutedNValues, int sizeDifference)
        {
            if (sizeDifference != 0)
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
        public List<List<string>> GenerateAllCombinationsOfTwoConditions(List<string> conditions)
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
        public void PrintSignificantProtein(List<ProteinRowInfo> allProteinInfo, double nValueThreshold,
            Dictionary<string, List<double>> allProteinObservedStatistics, string printSignificantProteinsLocation)
        {
            using (StreamWriter writetext = new StreamWriter(printSignificantProteinsLocation))
            {
                writetext.Write("ProteinID" + ", " + "PValue" + ", " + "LogFoldChange" + ", " + "Significance");
                writetext.WriteLine();
                foreach (string proteinID in allProteinObservedStatistics.Keys)
                {
                    if (allProteinObservedStatistics[proteinID][0] < nValueThreshold)
                    {
                        writetext.Write(proteinID + ", " + allProteinObservedStatistics[proteinID][1] + ", " 
                            + allProteinObservedStatistics[proteinID][2] + ", " + "Not Significant");
                        writetext.WriteLine();
                    }
                    else
                    {
                        writetext.Write(proteinID + ", " + allProteinObservedStatistics[proteinID][1] + ", "
                            + allProteinObservedStatistics[proteinID][2] + ", " + "Significant");
                        writetext.WriteLine();
                    }
                }

                //for (int i = 0; i < observedNValues.Count; i++)
                //{
                //    if (i == 0)
                //    {
                //        writetext.Write("ProteinID" + ", " + "PValue" + ", " + "LogFoldChange" + ", " + "Significance");
                //        writetext.WriteLine();
                //    }
                //    ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                //    if (observedNValues[i] < nValueThreshold)
                //    {
                //        writetext.Write(proteinRowInfo.ProteinID + ", " + observedPValues[i] + ", " + observedLogFoldChange[i]
                //            + ", " + "Not Significant");
                //        writetext.WriteLine();
                //    }
                //    else
                //    {
                //        writetext.Write(proteinRowInfo.ProteinID + ", " + observedPValues[i] + ", " + observedLogFoldChange[i]
                //            + ", " + "Significant");
                //        writetext.WriteLine();
                //    }
                //}
            }
        }

        public static void Main(string[] args)
        {
            RunProteinSignificanceClassifier proteinBasedSignificance = new RunProteinSignificanceClassifier();
            int numberSamples = 0;
            // Parse the ExperimentalDesign File to get info of samples and conditions they belong to
            Dictionary<string, List<string>> samplefileConditionRelation = proteinBasedSignificance.ExpermientalDesignParser("C:/Users/Anay/Desktop/UW Madison/Smith Lab/Spectra Data/ExperimentalDesign.tsv"
                , ref numberSamples);
            // get all conditions and pair them up for Significance classification
            List<string> allConditions = new List<string>(samplefileConditionRelation.Keys);
            List<List<string>> allTwoConditionCombinations = proteinBasedSignificance.GenerateAllCombinationsOfTwoConditions(allConditions);
            foreach (List<string> conditionPair in allTwoConditionCombinations)
            {
                string firstCondition = conditionPair[0];
                string secondCondition = conditionPair[1];
                double sOValue = 0.1;
                double meanFraction = 0.1;
                int maxSignificantCount = 0;

                while (meanFraction < 1)
                {
                    while (sOValue < 1)
                    {
                        for (int k = 0; k < numberSamples; k++)
                        {
                            proteinBasedSignificance = new RunProteinSignificanceClassifier();
                            //Declaring variables which will be generated after parsing QuantifiedPeptides file
                            List<ProteinRowInfo> allProteinInfo = new List<ProteinRowInfo>();
                            List<string> samplesFileNames = new List<string>();
                            int maxInvalidIntensityValues = k;
                            proteinBasedSignificance.ProteinDataParser(allProteinInfo, maxInvalidIntensityValues,
                                samplesFileNames, "C:/Users/Anay/Desktop/UW Madison/Smith Lab/Spectra Data/FlashLFQ_2020-04-26-17-39-35/QuantifiedProteins.tsv");

                            // imputes missing intensity values for each protein
                            ImputationProcess imputationProcess = new ImputationProcess();
                            imputationProcess.RunImputationProcess(allProteinInfo, samplesFileNames, meanFraction);

                            // Declaring variables which will be generated after T-Tests and Permutation Tests
                            List<double> observedNValues = new List<double>(); // will store observed N values
                            List<double> permutedNValues = new List<double>(); // will store permuted N values
                            StatisticalTests statisticalTests = new StatisticalTests();

                            // contains proteins and their observed N value, P Value and Log Fold Change 
                            Dictionary<string, List<double>> allProteinObservedStatistics = new Dictionary<string, List<double>>();
                            // Creating threads for Parallelizing code
                            ThreadPool.GetMaxThreads(out int workerThreadsCount, out int ioThreadsCount);
                            int[] threads = Enumerable.Range(0, workerThreadsCount).ToArray();
                            Parallel.ForEach(threads, (i) =>
                            {
                                // Compute observed and permuted N Values for each protein using T Tests and Permutation Testing
                                for (; i < allProteinInfo.Count; i += workerThreadsCount)
                                {
                                    ProteinRowInfo proteinRowInfo = allProteinInfo[i];
                                    Dictionary<string, double> samplesintensityData = proteinRowInfo.SamplesIntensityData;

                                    List<string> firstConditionAssociatedSamples = samplefileConditionRelation.GetValueOrDefault(firstCondition);
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

                                    // Compute observed N Values with the chosen pair of conditions using T-Tests and
                                    // store in observedNValues array
                                    List<double> proteinStatistics = statisticalTests.GetNValueUsingTTest(proteinFirstConditionIntensityValues, proteinSecondConditionIntensityValues,
                                        sOValue, false);

                                    // Compute permuted N Values with the chosen pair of conditions using T-Tests and 
                                    // store in permutedNValues array
                                    List<double> proteinPermutedNavlues = statisticalTests.GetNValueUsingPermutationtests(proteinFirstConditionIntensityValues, 
                                        proteinSecondConditionIntensityValues, sOValue);

                                    // add computed original and permuted statistics for the protein
                                    lock(allProteinObservedStatistics)
                                    {
                                        // add protein and its observed N value, P Value and Log Fold Change in that order
                                        allProteinObservedStatistics.Add(proteinRowInfo.ProteinID, new List<double>() { proteinStatistics[0],
                                        proteinStatistics[1], proteinStatistics[2]});
                                        observedNValues.Add(proteinStatistics[0]);
                                        foreach (double permutedNValue in proteinPermutedNavlues)
                                        {
                                            permutedNValues.Add(permutedNValue);
                                        }
                                    }
                                }

                            });

                            // makes the permuted N values list and the observed N Values list of the same size
                            proteinBasedSignificance.ResizePermutedArray(permutedNValues, permutedNValues.Count() - observedNValues.Count());

                            // get the threshold at which we will filter out the significant proteins
                            double nValueThreshold = statisticalTests.calculateNvaluethreshold(observedNValues, permutedNValues, 0.05);

                            // determine number of signifcant proteins detected
                            int newSignificantCount = observedNValues.Count(x => x >= nValueThreshold);
                            if (newSignificantCount > maxSignificantCount)
                            {
                                maxSignificantCount = newSignificantCount;
                                proteinBasedSignificance.PrintSignificantProtein(allProteinInfo, nValueThreshold, allProteinObservedStatistics, 
                                    "C:/Users/Anay/Desktop/UW Madison/Smith Lab/Project 1/ConsoleApp1/ProteinBasedSignificanceModified.csv");
                                Console.WriteLine("Sig Count - " + maxSignificantCount + "meanFraction - " + meanFraction + "sOValue - "
                                    + sOValue + "k - " + k);
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
}