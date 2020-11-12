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
        //ProteinID is never used. I'm wondering if this object isn't necessary and you could simply replace it with a Dictionary?
        public ProteinRowInfo()
        {
            SamplesIntensityData = new Dictionary<string, double>();
        }

        /// <summary>
        /// Stores detected intensity values of protein in each sample. 
        /// Maps condition's sample name to the intensity value
        /// </summary>
        public Dictionary<string, double> SamplesIntensityData { get; set; }

        /// <summary>
        /// Unique protein ID
        /// </summary>
        public string ProteinID { get; set; }
    }
}