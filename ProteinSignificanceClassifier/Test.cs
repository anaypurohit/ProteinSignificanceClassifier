using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Text;

namespace ProteinSignificanceClassifier
{
    [TestFixture]
    public class Test
    {
        /// <summary>
        /// These are the Unit Tests for the ImputationProcess class
        /// </summary>
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

        [TestCase]
        public void TestImputeData()
        {
            ProteinRowInfo expectedSampleprotein = new ProteinRowInfo();
            expectedSampleprotein.ProteinID = "A0A087WTI9";
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 24.607;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 23.763;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 20.783;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 20.165;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 20.639;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 17.223;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 21.564;
            expectedSampleprotein.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 20.469;


            List<string> sampleFileNames = new List<string>();
            sampleFileNames.Add("Intensity_B02_06_161103_A1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_07_161103_A2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_16_161103_A3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_17_161103_A4_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_24_161103_C1_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_09_161103_C2_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_14_161103_C3_HCD_OT_4ul");
            sampleFileNames.Add("Intensity_B02_19_161103_C4_HCD_OT_4ul");

            double[] sampleMissingCount = new double[8];
            sampleMissingCount[0] = 20;
            sampleMissingCount[1] = 7;
            sampleMissingCount[2] = 11;
            sampleMissingCount[3] = 0;
            sampleMissingCount[4] = 32;
            sampleMissingCount[5] = 2;
            sampleMissingCount[6] = 0;
            sampleMissingCount[7] = 5;

            int[] numberOfValidIntensityValuesInSample = new int[8];
            numberOfValidIntensityValuesInSample[0] = 120;
            numberOfValidIntensityValuesInSample[1] = 135;
            numberOfValidIntensityValuesInSample[2] = 100;
            numberOfValidIntensityValuesInSample[3] = 87;
            numberOfValidIntensityValuesInSample[4] = 189;
            numberOfValidIntensityValuesInSample[5] = 50;
            numberOfValidIntensityValuesInSample[6] = 79;
            numberOfValidIntensityValuesInSample[7] = 30;

            double[] samplesMeanIntensityValue = new double[8];
            samplesMeanIntensityValue[0] = 24.4;
            samplesMeanIntensityValue[1] = 23.5;
            samplesMeanIntensityValue[2] = 22;
            samplesMeanIntensityValue[3] = 19.7;
            samplesMeanIntensityValue[4] = 21.75;
            samplesMeanIntensityValue[5] = 18.2;
            samplesMeanIntensityValue[6] = 23.21;
            samplesMeanIntensityValue[7] = 21.09;

            double[] samplesStandardDeviation = new double[8];
            samplesStandardDeviation[0] = 0.2;
            samplesStandardDeviation[1] = 1.1;
            samplesStandardDeviation[2] = 0.02;
            samplesStandardDeviation[3] = 0.7;
            samplesStandardDeviation[4] = 1.7;
            samplesStandardDeviation[5] = 0.8;
            samplesStandardDeviation[6] = 2.4;
            samplesStandardDeviation[7] = 0.15;

            ProteinRowInfo outputSampleprotein = new ProteinRowInfo();
            outputSampleprotein.ProteinID = "A0A087WTI9";
            outputSampleprotein.SamplesIntensityData["Intensity_B02_06_161103_A1_HCD_OT_4ul"] = 0;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_07_161103_A2_HCD_OT_4ul"] = 23.763;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_16_161103_A3_HCD_OT_4ul"] = 20.783;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_17_161103_A4_HCD_OT_4ul"] = 20.165;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_24_161103_C1_HCD_OT_4ul"] = 0;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_09_161103_C2_HCD_OT_4ul"] = 0;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_14_161103_C3_HCD_OT_4ul"] = 21.564;
            outputSampleprotein.SamplesIntensityData["Intensity_B02_19_161103_C4_HCD_OT_4ul"] = 20.469;
            double meanFraction = 0.2;

            ImputationProcess imputationProcess = new ImputationProcess();
            imputationProcess.ImputeData(outputSampleprotein, samplesMeanIntensityValue, samplesStandardDeviation,
                sampleFileNames, sampleMissingCount, numberOfValidIntensityValuesInSample, meanFraction, true);

            Assert.That(expectedSampleprotein.SamplesIntensityData, 
                Is.EqualTo(outputSampleprotein.SamplesIntensityData).Within(0.001));
        }

        /// <summary>
        /// These are the Unit Tests for the StatisticalTests class
        /// </summary>
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

        [TestCase]
        public void TestGetNValueUsingTTest()
        {
            double expectedNValue = -9.276;
            double expectedPValue = 0.137;
            double expectedLogFoldChange = -0.087;

            List<double> proteinFirstConditionIntensityValues = new List<double>();
            List<double> proteinSecondConditionIntensityValues = new List<double>();
            proteinFirstConditionIntensityValues.Add(25.535);
            proteinFirstConditionIntensityValues.Add(25.482);
            proteinFirstConditionIntensityValues.Add(25.308);
            proteinFirstConditionIntensityValues.Add(25.373);

            proteinSecondConditionIntensityValues.Add(25.370);
            proteinSecondConditionIntensityValues.Add(25.368);
            proteinSecondConditionIntensityValues.Add(25.359);
            proteinSecondConditionIntensityValues.Add(25.251);

            double sOValue = 0.3;

            StatisticalTests statisticalTests = new StatisticalTests();
            List<double> proteinStatistics = statisticalTests.GetNValueUsingTTest(proteinFirstConditionIntensityValues, 
                proteinSecondConditionIntensityValues, sOValue, false);

            Assert.AreEqual(expectedNValue, proteinStatistics[0], 0.001);
            Assert.AreEqual(expectedPValue, proteinStatistics[1], 0.001);
            Assert.AreEqual(expectedLogFoldChange, proteinStatistics[2], 0.001);
        }

        [TestCase]
        public void TestGetNValueUsingPermutationTests()
        {
            List<double> expectedPermutedNValues = new List<double>();
            expectedPermutedNValues.Add(-13.021);
            expectedPermutedNValues.Add(-12.269);
            expectedPermutedNValues.Add(-8.298);
            expectedPermutedNValues.Add(-12.119);
            expectedPermutedNValues.Add(-8.270);
            expectedPermutedNValues.Add(-8.161);
            expectedPermutedNValues.Add(-17.649);
            expectedPermutedNValues.Add(-19.749);
            expectedPermutedNValues.Add(-24.623);
            expectedPermutedNValues.Add(-20.303);
            expectedPermutedNValues.Add(-23.756);
            expectedPermutedNValues.Add(-20.594);


            List<double> proteinFirstConditionIntensityValues = new List<double>();
            List<double> proteinSecondConditionIntensityValues = new List<double>();
            proteinFirstConditionIntensityValues.Add(25.535);
            proteinFirstConditionIntensityValues.Add(25.482);
            proteinFirstConditionIntensityValues.Add(25.308);
            proteinFirstConditionIntensityValues.Add(25.373);

            proteinSecondConditionIntensityValues.Add(25.370);
            proteinSecondConditionIntensityValues.Add(25.368);
            proteinSecondConditionIntensityValues.Add(25.359);
            proteinSecondConditionIntensityValues.Add(25.251);

            
            double sOValue = 0.3;

            StatisticalTests statisticalTests = new StatisticalTests();
            List<double> outputPermutedNValues = statisticalTests.GetNValueUsingPermutationtests(proteinFirstConditionIntensityValues, 
                proteinSecondConditionIntensityValues, sOValue);

            Assert.That(expectedPermutedNValues, Is.EqualTo(outputPermutedNValues).Within(0.001));
        }


        /// <summary>
        /// These are the Unit Tests for the RunProteinSignificanceClassifier class
        /// </summary>
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
                    Assert.That(sampleprotein1.SamplesIntensityData,
                        Is.EqualTo(outputProtein.SamplesIntensityData).Within(0.001));

                }
                else if (outputProtein.ProteinID.Equals(sampleprotein2.ProteinID))
                {
                    Assert.That(sampleprotein2.SamplesIntensityData,
                        Is.EqualTo(outputProtein.SamplesIntensityData).Within(0.001));
                }
                else if (outputProtein.ProteinID.Equals(sampleprotein3.ProteinID))
                {
                    Assert.That(sampleprotein3.SamplesIntensityData,
                        Is.EqualTo(outputProtein.SamplesIntensityData).Within(0.001));
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

            List<double> expectedResizedArray = new List<double>() { 0, 1, 2, 3, 5, 6, 7, 9, 10, 12, 13, 14, 16, 17, 18 };

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


    }
}

