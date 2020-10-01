using System;
using System.Collections.Generic;
using System.Text;

namespace ProteinSignificanceClassifier
{
    /// <summary>
    /// Store information pertaining each unique detected protein
    /// </summary>
    public class ProteinRowInfo
    {
        /// <summary>
        /// Stores detected intensity values of protein in each sample. 
        /// Maps condition's sample name to the intensity value
        /// </summary>
        private Dictionary<string, double> samplesintensityData = new Dictionary<string, double>();

        /// <summary>
        /// Unique protein ID
        /// </summary>
        private string proteinID;


        public Dictionary<string, double> getSamplesIntensityValues()
        {
            return this.samplesintensityData;
        }

        public void setSamplesIntensityValues(string intensityName, double value)
        {
            this.samplesintensityData.Add(intensityName, value);
        }

        public string getProteinID()
        {
            return this.proteinID;
        }

        public void setProteinID(string proteinID)
        {
            this.proteinID = proteinID;
        }
    }
}
