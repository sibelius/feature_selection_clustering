/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    AIR.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.PKIDiscretize;

/**
<!-- globalinfo-start -->
* Attribute Interdependence Redundancy Measure.For more information see:<br/>
* <br/>
* Wai-Ho Au, Keith C. C. Chan, Andrew K. C. Wong, Yang Wang (2005). Attribute Clustering for Grouping, Selection, and Classification of Gene Expression Data. IEEE/ACM Transactions on Computational Biology and Bioinformatics. 2(2):83-101.
* <p/>
<!-- globalinfo-end -->
 *
<!-- technical-bibtex-start -->
* BibTeX:
* <pre>
* &#64;article{Au2005,
*    author = {Wai-Ho Au and Keith C. C. Chan and Andrew K. C. Wong and Yang Wang},
*    journal = {IEEE/ACM Transactions on Computational Biology and Bioinformatics},
*    number = {2},
*    pages = {83-101},
*    title = {Attribute Clustering for Grouping, Selection, and Classification of Gene Expression Data},
*    volume = {2},
*    year = {2005},
*    ISSN = {1545-5963}
* }
* </pre>
* <p/>
<!-- technical-bibtex-end -->
 *
 *
<!-- options-start -->
* Valid options are: <p/>
* 
* <pre> -F &lt;discretizer specification&gt;
*  Specify the discretizer algorithm for continuous features (default: weka.filters.unsupervised.attribute.PKIDiscretize)</pre>
* 
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class AIR
	extends PreBuiltDistanceMatrix
	implements DiscreteDistance, Serializable, TechnicalInformationHandler {

    /** for serialization */
    private static final long serialVersionUID = -8802935013143715943L;

    /** filter used for discretization of the dataset */
    protected Filter discFilter = new PKIDiscretize();

    /** flag to indicate if its needed to save the dataset discretized */
    protected boolean saveDiscretedData = false;

    /** discretized dataset */
    protected Instances discreteData = null;

    /**
     * Generate the distance matrix of a dataset
     * @throws Exception
     */
    public void makeDistanceMatrix() throws Exception {

	int nInstances = m_Instances.numInstances();
	int nAttributes = m_Instances.numAttributes();

	int classIndex = m_Instances.classIndex();

	if(classIndex >= 0){
	    m_DistMatrix = new double[nAttributes-1][nAttributes-1];
	}else{
	    m_DistMatrix = new double[nAttributes][nAttributes];
	}



	discFilter.setInputFormat(m_Instances);
	Instances cInstances = Filter.useFilter(m_Instances, discFilter);

	if (saveDiscreteData()) {
	    discreteData = cInstances;
	}

	for (int i = 0; i < nAttributes; ++i) {

	    //if the class is in the middle of the features, skip it
	    if(i == classIndex){
		continue;
	    }
	    
	    int ni = cInstances.attribute(i).numValues();
	    for (int j = i + 1; j < nAttributes; ++j) {
		//if the class is in the middle of the features, skip it
		if(j == classIndex){
		    continue;
		}

		int nj = cInstances.attribute(j).numValues();
		//contingency table
		double[][] counts = ContingencyTableExtended.makeContingencyTable(cInstances, i, j);
		double[] sumi = new double[ni];
		double[] sumj = new double[nj];

		int ii, jj;
		// values sums
		for (ii = 0; ii < ni; ii++) {
		    sumi[ii] = 0.0;

		    for (jj = 0; jj < nj; jj++) {
			sumi[ii] += counts[ii][jj];
		    }

		}

		for (jj = 0; jj < nj; jj++) {
		    sumj[jj] = 0.0;

		    for (ii = 0; ii < ni; ii++) {
			sumj[jj] += counts[ii][jj];
		    }
		}

		//assumes no missing values
		double mutualInformation = 0;
		double jointEntropy = 0;
		double probCond = 0;
		for (ii = 0; ii < ni; ++ii) {
		    for (jj = 0; jj < nj; ++jj) {
			probCond = counts[ii][jj] / nInstances;
			if (probCond == 0) {
			    continue;
			}
			mutualInformation += probCond * Utils.log2(probCond / ((sumi[ii] / nInstances) * (sumj[jj] / nInstances)));
			jointEntropy -= probCond * Utils.log2(probCond);

		    }
		}

		//this makes the features that are after the class reduce their indexes by one as done by the invertDataset filter
		int idxI = i;
		int idxJ = j;
		if(idxI > classIndex){
		    --idxI;
		}
		if(idxJ > classIndex){
		    --idxJ;
		}
		m_DistMatrix[idxI][idxJ] =
			m_DistMatrix[idxJ][idxI] = 1 - (mutualInformation / jointEntropy);
		//System.err.println("["+idxI+"]["+idxJ+"] "+m_DistMatrix[idxI][idxJ]);
	    }
	}
    }

    /**
     * Return the distance between two features
     * @param first a feature inverted using the InvertDataset filter
     * @param second a feature inverted using the InvertDataset filter
     * @return the distance between the features
     */
    public double distance(Instance first, Instance second) {
	return m_DistMatrix[(int) first.value(0)][(int) second.value(0)];
    }


    /**
     * Returns an enumeration describing the available options.
     * @return an enumeration of all the available options.
     *
     **/
    public Enumeration<Option> listOptions() {
         Vector<Option> newVector = new Vector<Option>();
         newVector.addElement(new Option("\tSpecify the discretizer algorithm for continuous " +
		"features (default: weka.filters.unsupervised.attribute.PKIDiscretize)",
                "F", 1, "-F <discretizer specification>"));
        return newVector.elements();
    }

    /**
     * Parses a given list of options.
     *
    <!-- options-start -->
    * Valid options are: <p/>
    * 
    * <pre> -F &lt;discretizer specification&gt;
    *  Specify the discretizer algorithm for continuous features (default: weka.filters.unsupervised.attribute.PKIDiscretize)</pre>
    * 
    <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     *
     **/
    public void setOptions(String[] options) throws Exception {
        String optionString = Utils.getOption('F', options);
	String [] discretizerSpec = Utils.splitOptions(optionString);
        if (optionString.length() != 0) {
	    String discretizerName = discretizerSpec[0];
	    discretizerSpec[0] = "";
            setDiscretizer((Filter) Utils.forName(Filter.class, discretizerName, discretizerSpec));
        }
    }

    /**
     * Gets the current settings.
     * @return an array of strings suitable for passing to setOptions()
     */
    public String[] getOptions() {
	Vector<String> options = new Vector<String>();

        options.add("-F");
	if( discFilter instanceof OptionHandler){
	    OptionHandler df = (OptionHandler)discFilter;
	    options.add(discFilter.getClass().getName()+ " " + Utils.joinOptions(df.getOptions()));
	}else{
	    options.add(discFilter.getClass().getName());
	}
        return options.toArray(new String[options.size()]);
    }

    /**
     * Reset options to default values
     */
    private void resetOptions() {
	discFilter = new PKIDiscretize();
    }

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "Attribute Interdependence Redundancy Measure."
		+ "For more information see:\n\n" + getTechnicalInformation().toString();
    }

    /**
     * Returns an instance of a TechnicalInformation object, containing
     * detailed information about the technical background of this class,
     * e.g., paper reference or book this class is based on.
     *
     * @return the technical information about this class
     */
    public TechnicalInformation getTechnicalInformation() {
	TechnicalInformation result;

        result = new TechnicalInformation(Type.ARTICLE);
        result.setValue(Field.AUTHOR, "Wai-Ho Au and Keith C. C. Chan and Andrew K. C. Wong and Yang Wang");
        result.setValue(Field.YEAR, "2005");
        result.setValue(Field.TITLE, "Attribute Clustering for Grouping, Selection, and " +
		"Classification of Gene Expression Data");
        result.setValue(Field.JOURNAL, "IEEE/ACM Transactions on Computational Biology and Bioinformatics");
        result.setValue(Field.VOLUME, "2");
	result.setValue(Field.NUMBER, "2");
	result.setValue(Field.ISSN, "1545-5963");
        result.setValue(Field.PAGES, "83-101");

	return result;
    }

    /**
     * Sets the filter used to discretize the dataset
     * @param discretizer discretize filter
     */
    public void setDiscretizer(Filter discretizer) {
	discFilter = discretizer;
    }

    /**
     * Get the filter used to discretize the dataset
     * @return discretize filter
     */
    public Filter getDiscretizer() {
	return discFilter;
    }

    /**
     * Tip text for this property
     * @return tip text for the GUI
     */
    public String discretizerTipText(){
	return "Filter used to discretize the dataset";
    }
    /**
     * Get the dataset discretized
     * @return the discretized dataset
     */
    public Instances getDiscreteData() {
	return discreteData;
    }

    /**
     * Set the flag for saving the discrete dataset
     * @param save flag updated
     */
    public void saveDiscreteData(boolean save) {
	saveDiscretedData = save;
    }

    /**
     * Get the flag which indicates if the discrete dataset is saved
     * @return true if discrete data is saved
     */
    public boolean saveDiscreteData() {
	return saveDiscretedData;
    }
}

