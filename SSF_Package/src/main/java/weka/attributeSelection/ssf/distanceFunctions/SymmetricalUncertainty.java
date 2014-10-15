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
 *    SymmetricalUncertainty.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

import weka.core.ContingencyTables;
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
* Symmetrical Uncertainty. For more informatio see: William H. Press, Brian P. Flannery, Saul A. Teukolsky, William T. Vetterling (1990). Numerical Recipes in C: The Art of Scientific Computing. .
* <p/>
<!-- globalinfo-end -->
 *
<!-- technical-bibtex-start -->
* BibTeX:
* <pre>
* &#64;book{Press1990,
*    author = {William H. Press and Brian P. Flannery and Saul A. Teukolsky and William T. Vetterling},
*    title = {Numerical Recipes in C: The Art of Scientific Computing},
*    year = {1990}
* }
* </pre>
* <p/>
<!-- technical-bibtex-end -->
 *
 *
<!-- options-start -->
* Valid options are: <p/>
* 
* <pre> -F &lt;discretizer spec&gt;
*  Specify the discretizer algorithm for continuous features</pre>
* 
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class SymmetricalUncertainty
        extends PreBuiltDistanceMatrix
        implements DiscreteDistance, Serializable, TechnicalInformationHandler {

    /** for serialization */
    private static final long serialVersionUID = 3175917369991542434L;

    /** filter used for discretization */
    protected Filter discFilter = new PKIDiscretize();

    /** flag to indicate if the discrete dataset should be saved */
    protected boolean saveDiscretedData;

    /** saved discrete dataset */
    protected Instances discreteData = null;

    /**
     * Generate the distance matrix of a dataset
     * @throws Exception
     */
    @Override
    protected void makeDistanceMatrix() throws Exception {
        //avoid remaking calculations
        if (m_DistMatrix != null) {
            return;
        }
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
            
            for (int j = i + 1; j < nAttributes; ++j) {

		//if the class is in the middle of the features, skip it
		if(j == classIndex){
		    continue;
		}
                //contingency table
                double[][] counts = ContingencyTableExtended.makeContingencyTable(cInstances, i, j);

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
                        m_DistMatrix[idxJ][idxI] = 1 - ContingencyTables.symmetricalUncertainty(counts);
                //System.out.println("["+i+"]["+j+"]"+m_DistMatrix[i][j]);
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
     * Returns an enumeration describing the available options.
     * @return an enumeration of all the available options.
     *
     **/
    public Enumeration<Option> listOptions() {
         Vector<Option> newVector = new Vector<Option>();

        newVector.addElement(new Option("\tSpecify the discretizer algorithm for continuous features",
                "F", 1, "-F <discretizer spec>"));
        return newVector.elements();
    }

    /**
     * Parses a given list of options.
     *
    <!-- options-start -->
    * Valid options are: <p/>
    * 
    * <pre> -F &lt;discretizer spec&gt;
    *  Specify the discretizer algorithm for continuous features</pre>
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

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "Symmetrical Uncertainty. For more informatio see: " +
		getTechnicalInformation().toString();
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

        result = new TechnicalInformation(Type.BOOK);
        result.setValue(Field.AUTHOR, "William H. Press and Brian P. Flannery and " +
		"Saul A. Teukolsky and William T. Vetterling");
        result.setValue(Field.YEAR, "1990");
        result.setValue(Field.TITLE, "Numerical Recipes in C: The Art of Scientific Computing");

	return result;
    }

}

