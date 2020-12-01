using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Text;

namespace ProteinSignificanceClassifier
{
    [TestFixture]
    public class Test
    {
        [TestCase]
        public void TestExperimentalDesignParserConditionSampleMapping()
        {
            RunProteinSignificanceClassifier runProteinSignificanceClassifier = new RunProteinSignificanceClassifier();

            Dictionary<string, List<string>> expectedSamplefileConditionRelation = new Dictionary<string, List<string>>();
            expectedSamplefileConditionRelation.Add("Ecoli 1x", new List<string>() { "Intensity_B02_06_161103_A1_HCD_OT_4ul",
            "Intensity_B02_07_161103_A2_HCD_OT_4ul", "Intensity_B02_16_161103_A3_HCD_OT_4ul", "Intensity_B02_17_161103_A4_HCD_OT_4ul"});
            expectedSamplefileConditionRelation.Add("Ecoli 2x", new List<string>() { "Intensity_B02_09_161103_C2_HCD_OT_4ul",
            "Intensity_B02_14_161103_C3_HCD_OT_4ul", "Intensity_B02_19_161103_C4_HCD_OT_4ul", "Intensity_B02_24_161103_C1_HCD_OT_4ul"});

            int resultNumberSamples = 0;
            Dictionary<string, List<string>> resultSamplefileConditionRelation =
                runProteinSignificanceClassifier.ExpermientalDesignParser("TestExperimentalDesign.tsv", ref resultNumberSamples);

            CollectionAssert.AreEqual(expectedSamplefileConditionRelation, resultSamplefileConditionRelation);
        }

        [TestCase]
        public void TestExperimentalDesignParserSampleNumbers()
        {
            RunProteinSignificanceClassifier runProteinSignificanceClassifier = new RunProteinSignificanceClassifier();

            int expectedNumberSamples = 8;

            int resultNumberSamples = 0;
            Dictionary<string, List<string>> resultSamplefileConditionRelation =
                runProteinSignificanceClassifier.ExpermientalDesignParser("TestExperimentalDesign.tsv", ref resultNumberSamples);

            Assert.AreEqual(expectedNumberSamples, resultNumberSamples);
        }

        [TestCase]
        public void TestProteinDataParserSampleFileNames()
        {
            List<string> expectedSampleFileNames = new List<string>();
            expectedSampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            expectedSampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            int maxAllowedMissingValues = 2;
            List<ProteinRowInfo> outputAllProteinInfo = new List<ProteinRowInfo>();
            List<string> outputSamplesFileNames = new List<string>();
            RunProteinSignificanceClassifier runProteinSignificanceClassifier = new RunProteinSignificanceClassifier();
            runProteinSignificanceClassifier.ProteinDataParser(outputAllProteinInfo, maxAllowedMissingValues,
                outputSamplesFileNames, "TestQuantifiedProteins.tsv");

            CollectionAssert.AreEqual(expectedSampleFileNames, outputSamplesFileNames);
        }

        [TestCase]
        public void TestProteinDataParserProteinIntensityMapping()
        {
            List<string> sampleFileNames = new List<string>();
            sampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            List<ProteinRowInfo> expectedAllProteinInfo = new List<ProteinRowInfo>();
            ProteinRowInfo sampleprotein1 = new ProteinRowInfo();
            sampleprotein1.ProteinID = "A0A023KES1";
            sampleprotein1.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 22.192;
            sampleprotein1.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 22.729;
            sampleprotein1.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 22.347;
            sampleprotein1.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 22.397;
            sampleprotein1.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 23.489;
            sampleprotein1.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 23.180;
            sampleprotein1.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 23.293;
            sampleprotein1.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 23.271;

            ProteinRowInfo sampleprotein2 = new ProteinRowInfo();
            sampleprotein2.ProteinID = "A0A024R4E5";
            sampleprotein2.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 25.535;
            sampleprotein2.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 25.482;
            sampleprotein2.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 25.308;
            sampleprotein2.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 25.373;
            sampleprotein2.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 25.370;
            sampleprotein2.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 25.368;
            sampleprotein2.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 25.359;
            sampleprotein2.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 25.251;

            ProteinRowInfo sampleprotein3 = new ProteinRowInfo();
            sampleprotein3.ProteinID = "A0A087WTI9";
            sampleprotein3.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein3.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 23.763;
            sampleprotein3.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 20.783;
            sampleprotein3.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 20.165;
            sampleprotein3.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 21.7858;
            sampleprotein3.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 0;
            sampleprotein3.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 21.564;
            sampleprotein3.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 20.469;

            expectedAllProteinInfo.Add(sampleprotein1);
            expectedAllProteinInfo.Add(sampleprotein2);
            expectedAllProteinInfo.Add(sampleprotein3);


            int maxAllowedMissingValues = 2;
            List<ProteinRowInfo> outputAllProteinInfo = new List<ProteinRowInfo>();
            List<string> outputSamplesFileNames = new List<string>();
            RunProteinSignificanceClassifier runProteinSignificanceClassifier = new RunProteinSignificanceClassifier();
            runProteinSignificanceClassifier.ProteinDataParser(outputAllProteinInfo, maxAllowedMissingValues,
                outputSamplesFileNames, "TestQuantifiedProteins.tsv");

            foreach (ProteinRowInfo outputProtein in outputAllProteinInfo)
            {
                if (outputProtein.ProteinID.Equals(sampleprotein1.ProteinID))
                {
                    foreach (string sampleName in sampleFileNames)
                    {
                        Assert.AreEqual(outputProtein.SamplesIntensityData[sampleName],
                            sampleprotein1.SamplesIntensityData[sampleName], 0.001);
                    }
                }
                else if (outputProtein.ProteinID.Equals(sampleprotein2.ProteinID))
                {
                    foreach (string sampleName in outputSamplesFileNames)
                    {
                        Assert.AreEqual(outputProtein.SamplesIntensityData[sampleName],
                            sampleprotein2.SamplesIntensityData[sampleName], 0.001);
                    }
                }
                else if (outputProtein.ProteinID.Equals(sampleprotein3.ProteinID))
                {
                    foreach (string sampleName in outputSamplesFileNames)
                    {
                        Assert.AreEqual(outputProtein.SamplesIntensityData[sampleName],
                            sampleprotein3.SamplesIntensityData[sampleName], 0.001);
                    }
                }
                else
                {
                    Assert.Fail();
                }
            }
        }

        [TestCase]
        public void TestResizePermutedArray()
        {

            List<double> expectedResizedArray = new List<double>() { 0, 1, 2, 3, 5, 6, 7, 9, 10, 12, 13, 14, 16, 17, 18};

            List<double> outputResizedArray = new List<double>();
            for (int i = 0; i < 20; i++)
            {
                outputResizedArray.Add(i);
            }
            int sizeDifference = 5;
            RunProteinSignificanceClassifier runProteinSignificanceClassifier = new RunProteinSignificanceClassifier();
            runProteinSignificanceClassifier.ResizePermutedArray(outputResizedArray, sizeDifference);

            CollectionAssert.AreEqual(expectedResizedArray, outputResizedArray);
        }

        [TestCase]
        public void TestRPCGenerateAllCombinationsOfTwoConditions()
        {
            RunProteinSignificanceClassifier runProteinSignificanceClassifier = new RunProteinSignificanceClassifier();
            List<string> testConditions = new List<string>();
            testConditions.Add("condition1");
            testConditions.Add("condition2");
            testConditions.Add("condition3");
            testConditions.Add("condition4");

            List<List<string>> expectedConditionPairs = new List<List<string>>();
            expectedConditionPairs.Add(new List<string>() { "condition1", "condition2" });
            expectedConditionPairs.Add(new List<string>() { "condition1", "condition3" });
            expectedConditionPairs.Add(new List<string>() { "condition1", "condition4" });
            expectedConditionPairs.Add(new List<string>() { "condition2", "condition3" });
            expectedConditionPairs.Add(new List<string>() { "condition2", "condition4" });
            expectedConditionPairs.Add(new List<string>() { "condition3", "condition4" });

            List<List<string>> resultConditionPairs = runProteinSignificanceClassifier.GenerateAllCombinationsOfTwoConditions(testConditions);
            CollectionAssert.AreEqual(expectedConditionPairs, resultConditionPairs);
        }

        // Imputation Process Tests

        [TestCase]
        public void TestCalculateMissingIntensityValuesInSample()
        {
            double[] expectedSampleMissingCount = new double[8];
            expectedSampleMissingCount[0] = 1;
            expectedSampleMissingCount[1] = 0;
            expectedSampleMissingCount[2] = 0;
            expectedSampleMissingCount[3] = 0;
            expectedSampleMissingCount[4] = 1;
            expectedSampleMissingCount[5] = 0;
            expectedSampleMissingCount[6] = 0;
            expectedSampleMissingCount[7] = 1;

            List<string> sampleFileNames = new List<string>();
            sampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            ProteinRowInfo sampleprotein1 = new ProteinRowInfo();
            sampleprotein1.ProteinID = "A0A023KES1";
            sampleprotein1.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 22.729;
            sampleprotein1.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 22.347;
            sampleprotein1.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 22.397;
            sampleprotein1.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 23.180;
            sampleprotein1.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 23.293;
            sampleprotein1.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;


            ProteinRowInfo sampleprotein2 = new ProteinRowInfo();
            sampleprotein2.ProteinID = "A0A024R4E5";
            sampleprotein2.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 25.535;
            sampleprotein2.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 25.370;
            sampleprotein2.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 24;
            sampleprotein2.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 25.359;
            sampleprotein2.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;


            double[] outputSampleMissingCount = new double[8];
            ImputationProcess imputationProcess = new ImputationProcess();
            imputationProcess.CalculateMissingIntensityValuesInSample(sampleprotein1,
                outputSampleMissingCount, sampleFileNames);

            CollectionAssert.AreEqual(expectedSampleMissingCount, outputSampleMissingCount);

            expectedSampleMissingCount[0] = 2;
            expectedSampleMissingCount[1] = 0;
            expectedSampleMissingCount[2] = 1;
            expectedSampleMissingCount[3] = 1;
            expectedSampleMissingCount[4] = 1;
            expectedSampleMissingCount[5] = 0;
            expectedSampleMissingCount[6] = 0;
            expectedSampleMissingCount[7] = 2;

            imputationProcess.CalculateMissingIntensityValuesInSample(sampleprotein2,
            outputSampleMissingCount, sampleFileNames);

            CollectionAssert.AreEqual(expectedSampleMissingCount, outputSampleMissingCount);
        }

        [TestCase]
        public void TestCalculateNumberOfIntensityValuesOfSample()
        {
            int[] expectedSampleIntensityValuesCount = new int[8];
            expectedSampleIntensityValuesCount[0] = 0;
            expectedSampleIntensityValuesCount[1] = 1;
            expectedSampleIntensityValuesCount[2] = 1;
            expectedSampleIntensityValuesCount[3] = 1;
            expectedSampleIntensityValuesCount[4] = 0;
            expectedSampleIntensityValuesCount[5] = 1;
            expectedSampleIntensityValuesCount[6] = 1;
            expectedSampleIntensityValuesCount[7] = 0;

            List<string> sampleFileNames = new List<string>();
            sampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            ProteinRowInfo sampleprotein1 = new ProteinRowInfo();
            sampleprotein1.ProteinID = "A0A023KES1";
            sampleprotein1.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 22.729;
            sampleprotein1.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 22.347;
            sampleprotein1.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 22.397;
            sampleprotein1.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 23.180;
            sampleprotein1.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 23.293;
            sampleprotein1.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;


            ProteinRowInfo sampleprotein2 = new ProteinRowInfo();
            sampleprotein2.ProteinID = "A0A024R4E5";
            sampleprotein2.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 25.535;
            sampleprotein2.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 25.370;
            sampleprotein2.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 24;
            sampleprotein2.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 25.359;
            sampleprotein2.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;

            int[] outputSampleIntensityValuesCount = new int[8];
            double[] outputSampleAllIntensityValuesSum = new double[8];
            ImputationProcess imputationProcess = new ImputationProcess();
            imputationProcess.CalculateNumberAndSumOfIntensityValuesOfSample(sampleprotein1,
                outputSampleAllIntensityValuesSum, sampleFileNames, outputSampleIntensityValuesCount);

            CollectionAssert.AreEqual(expectedSampleIntensityValuesCount, outputSampleIntensityValuesCount);

            expectedSampleIntensityValuesCount[0] = 0;
            expectedSampleIntensityValuesCount[1] = 2;
            expectedSampleIntensityValuesCount[2] = 1;
            expectedSampleIntensityValuesCount[3] = 1;
            expectedSampleIntensityValuesCount[4] = 1;
            expectedSampleIntensityValuesCount[5] = 2;
            expectedSampleIntensityValuesCount[6] = 2;
            expectedSampleIntensityValuesCount[7] = 0;

            imputationProcess.CalculateNumberAndSumOfIntensityValuesOfSample(sampleprotein2,
                outputSampleAllIntensityValuesSum, sampleFileNames, outputSampleIntensityValuesCount);

            CollectionAssert.AreEqual(expectedSampleIntensityValuesCount, outputSampleIntensityValuesCount);
        }

        [TestCase]
        public void TestCalculateSumOfIntensityValuesOfSample()
        {
            double[] expectedSampleAllIntensityValuesSum = new double[8];
            expectedSampleAllIntensityValuesSum[0] = 0;
            expectedSampleAllIntensityValuesSum[1] = 22.729;
            expectedSampleAllIntensityValuesSum[2] = 22.347;
            expectedSampleAllIntensityValuesSum[3] = 22.397;
            expectedSampleAllIntensityValuesSum[4] = 0;
            expectedSampleAllIntensityValuesSum[5] = 23.180;
            expectedSampleAllIntensityValuesSum[6] = 23.293;
            expectedSampleAllIntensityValuesSum[7] = 0;

            List<string> sampleFileNames = new List<string>();
            sampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            ProteinRowInfo sampleprotein1 = new ProteinRowInfo();
            sampleprotein1.ProteinID = "A0A023KES1";
            sampleprotein1.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 22.729;
            sampleprotein1.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 22.347;
            sampleprotein1.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 22.397;
            sampleprotein1.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 23.180;
            sampleprotein1.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 23.293;
            sampleprotein1.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;

            ProteinRowInfo sampleprotein2 = new ProteinRowInfo();
            sampleprotein2.ProteinID = "A0A024R4E5";
            sampleprotein2.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 25.535;
            sampleprotein2.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 25.370;
            sampleprotein2.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 24;
            sampleprotein2.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 25.359;
            sampleprotein2.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;

            int[] outputSampleIntensityValuesCount = new int[8];
            double[] outputSampleAllIntensityValuesSum = new double[8];
            ImputationProcess imputationProcess = new ImputationProcess();
            imputationProcess.CalculateNumberAndSumOfIntensityValuesOfSample(sampleprotein1,
                outputSampleAllIntensityValuesSum, sampleFileNames, outputSampleIntensityValuesCount);
            CollectionAssert.AreEqual(expectedSampleAllIntensityValuesSum, outputSampleAllIntensityValuesSum);

            expectedSampleAllIntensityValuesSum[0] = 0;
            expectedSampleAllIntensityValuesSum[1] = 22.729 + 25.535;
            expectedSampleAllIntensityValuesSum[2] = 22.347;
            expectedSampleAllIntensityValuesSum[3] = 22.397;
            expectedSampleAllIntensityValuesSum[4] = 25.370;
            expectedSampleAllIntensityValuesSum[5] = 23.180 + 24;
            expectedSampleAllIntensityValuesSum[6] = 23.293 + 25.359;
            expectedSampleAllIntensityValuesSum[7] = 0;

            imputationProcess.CalculateNumberAndSumOfIntensityValuesOfSample(sampleprotein2,
                outputSampleAllIntensityValuesSum, sampleFileNames, outputSampleIntensityValuesCount);
            CollectionAssert.AreEqual(expectedSampleAllIntensityValuesSum, outputSampleAllIntensityValuesSum);
        }

        [TestCase]
        public void TestCalculateSampleStandardDeviationNumerator()
        {
            List<string> sampleFileNames = new List<string>();
            sampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            ProteinRowInfo sampleprotein1 = new ProteinRowInfo();
            sampleprotein1.ProteinID = "A0A023KES1";
            sampleprotein1.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 22.729;
            sampleprotein1.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 22.347;
            sampleprotein1.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 22.397;
            sampleprotein1.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 0;
            sampleprotein1.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 23.180;
            sampleprotein1.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 23.293;
            sampleprotein1.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;

            ProteinRowInfo sampleprotein2 = new ProteinRowInfo();
            sampleprotein2.ProteinID = "A0A024R4E5";
            sampleprotein2.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 25.535;
            sampleprotein2.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 0;
            sampleprotein2.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 25.370;
            sampleprotein2.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 24;
            sampleprotein2.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 25.359;
            sampleprotein2.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 0;

            double[] samplesMeanIntensityValue = new double[8];
            samplesMeanIntensityValue[0] = 0;
            samplesMeanIntensityValue[1] = (22.729 + 25.535)/2;
            samplesMeanIntensityValue[2] = 22.347/2;
            samplesMeanIntensityValue[3] = 22.397/2;
            samplesMeanIntensityValue[4] = 25.370/2;
            samplesMeanIntensityValue[5] = (23.180 + 24)/2;
            samplesMeanIntensityValue[6] = (23.293 + 25.359)/2;
            samplesMeanIntensityValue[7] = 0;

            double[] expectedSamplesStandardDeviationNumerators = new double[8];
            expectedSamplesStandardDeviationNumerators[0] = 0;
            expectedSamplesStandardDeviationNumerators[1] = Math.Pow(22.729 - samplesMeanIntensityValue[1], 2);
            expectedSamplesStandardDeviationNumerators[2] = Math.Pow(22.347 - samplesMeanIntensityValue[2], 2);
            expectedSamplesStandardDeviationNumerators[3] = Math.Pow(22.397 - samplesMeanIntensityValue[3], 2);
            expectedSamplesStandardDeviationNumerators[4] = 0;
            expectedSamplesStandardDeviationNumerators[5] = Math.Pow(23.180 - samplesMeanIntensityValue[5], 2);
            expectedSamplesStandardDeviationNumerators[6] = Math.Pow(23.293 - samplesMeanIntensityValue[6], 2);
            expectedSamplesStandardDeviationNumerators[7] = 0;

            double[] outputSamplesStandardDeviationNumerators = new double[8];
            ImputationProcess imputationProcess = new ImputationProcess();
            imputationProcess.CalculateSampleStandardDeviationNumerator(sampleprotein1, outputSamplesStandardDeviationNumerators,
                sampleFileNames, samplesMeanIntensityValue);
            CollectionAssert.AreEqual(expectedSamplesStandardDeviationNumerators, outputSamplesStandardDeviationNumerators);

            expectedSamplesStandardDeviationNumerators[0] = 0;
            expectedSamplesStandardDeviationNumerators[1] = Math.Pow(22.729 - samplesMeanIntensityValue[1], 2)
                + Math.Pow(25.535 - samplesMeanIntensityValue[1], 2);
            expectedSamplesStandardDeviationNumerators[2] = Math.Pow(22.347 - samplesMeanIntensityValue[2], 2);
            expectedSamplesStandardDeviationNumerators[3] = Math.Pow(22.397 - samplesMeanIntensityValue[3], 2);
            expectedSamplesStandardDeviationNumerators[4] = Math.Pow(25.370 - samplesMeanIntensityValue[4], 2);
            expectedSamplesStandardDeviationNumerators[5] = Math.Pow(23.180 - samplesMeanIntensityValue[5], 2)
                + Math.Pow(24 - samplesMeanIntensityValue[5], 2);
            expectedSamplesStandardDeviationNumerators[6] = Math.Pow(23.293 - samplesMeanIntensityValue[6], 2)
                + Math.Pow(25.359 - samplesMeanIntensityValue[6], 2);
            expectedSamplesStandardDeviationNumerators[7] = 0;

            imputationProcess.CalculateSampleStandardDeviationNumerator(sampleprotein2, outputSamplesStandardDeviationNumerators,
                sampleFileNames, samplesMeanIntensityValue);
            CollectionAssert.AreEqual(expectedSamplesStandardDeviationNumerators, outputSamplesStandardDeviationNumerators);
        }

        // still have to test the impute method in order to actually impute values

        // statistical tests

        [TestCase]
        public void TestStatisticalTestsGenerateAllCombinationsOfTwoConditions()
        {
            StatisticalTests statisticalTests = new StatisticalTests();

            List<double> testIndices = new List<double>();
            testIndices.Add(1);
            testIndices.Add(2);
            testIndices.Add(3);
            testIndices.Add(4);

            List<List<int>> expectedIndicesPairs = new List<List<int>>();
            expectedIndicesPairs.Add(new List<int>() { 0, 1 });
            expectedIndicesPairs.Add(new List<int>() { 0, 2 });
            expectedIndicesPairs.Add(new List<int>() { 0, 3 });
            expectedIndicesPairs.Add(new List<int>() { 1, 2 });
            expectedIndicesPairs.Add(new List<int>() { 1, 3 });
            expectedIndicesPairs.Add(new List<int>() { 2, 3 });

            List<List<int>> resultIndicesPairs = statisticalTests.GenerateAllCombinationsOfTwoIndices(testIndices);
            CollectionAssert.AreEqual(expectedIndicesPairs, resultIndicesPairs);
        }

        [TestCase]
        public void TestCalculateProteinIntensityValuesStandardDeviation()
        {
            List<double> testIntensityValues = new List<double>();
            testIntensityValues.Add(24.75);
            testIntensityValues.Add(25.15);
            testIntensityValues.Add(28.35);
            testIntensityValues.Add(21.95);
            double intensityValuesMean = 25.05;
            double expectedStandardDev = 2.269;

            StatisticalTests statisticalTests = new StatisticalTests();
            double outputStandardDev = statisticalTests.CalculateProteinIntensityValuesStandardDeviation
                (testIntensityValues, intensityValuesMean);

            Assert.AreEqual(expectedStandardDev, outputStandardDev, 0.001);
        }

        [TestCase]
        public void TestCalculateNvaluethreshold()
        {

            double expectedNValueThreshold = 74;

            List<double> observedNValues = new List<double>();
            for (int i = 0; i < 100; i++)
            {
                observedNValues.Add(i + 100);
            }

            List<double> permutedNValues = new List<double>();
            for (int i = 0; i < 100; i++)
            {
                permutedNValues.Add(i);
            }

            double FDR = 0.25;
            StatisticalTests statisticalTests = new StatisticalTests();
            double outputNValueThrehold = statisticalTests.calculateNvaluethreshold(observedNValues, 
                permutedNValues, FDR);

            Assert.AreEqual(expectedNValueThreshold, outputNValueThrehold, 0.001);
        }
    }
}

