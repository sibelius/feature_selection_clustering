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
 *    S3F.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection;

import java.io.Serializable;
import java.util.BitSet;
import java.util.Enumeration;
import java.util.Vector;
import weka.attributeSelection.ssf.clusterers.KMedoids;
import weka.attributeSelection.ssf.clusterers.KMedoidsSupervised;
import weka.attributeSelection.ssf.distanceFunctions.DistanceFunctionForSSF;
import weka.attributeSelection.ssf.distanceFunctions.SupervisedDistanceFunctionForSSF;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.SelectedTag;
import weka.core.Tag;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;

/**
<!-- globalinfo-start -->
* S3F:<br/>
* <br/>
* Select features from an induced partition of the feature set considering the correlation between features and class. See:<br/>
* <br/>
* Thiago F. Covões, Eduardo R. Hruschka (2011). Towards Improving Cluster-Based Feature Selection with a Simplified Silhouette Filter. Information Sciences. 181(18):3766-3782.
* <p/>
<!-- globalinfo-end -->
 *
<!-- technical-bibtex-start -->
* BibTeX:
* <pre>
* &#64;article{Covões2011,
*    author = {Thiago F. Covões and Eduardo R. Hruschka},
*    journal = {Information Sciences},
*    number = {18},
*    pages = {3766-3782},
*    title = {Towards Improving Cluster-Based Feature Selection with a Simplified Silhouette Filter},
*    volume = {181},
*    year = {2011},
*    ISSN = {0020-0255}
* }
* </pre>
* <p/>
<!-- technical-bibtex-end -->
 *
 *
<!-- options-start -->
* Valid options are: <p/>
* 
* <pre> -S &lt;num&gt;
*  Specify the strategy to select features from induced clusters (0 for medoids, and 1 for medoids and frontier).</pre>
* 
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class S3F
	extends ASEvaluation
	implements SSFInterface, OptionHandler, Serializable, TechnicalInformationHandler {

    /** for serialization */
    private static final long serialVersionUID = -1181782685346717350L;

    /** m_Strategy that will be used for picking attributes from the attribute clusters */
    protected int m_Strategy = STRATEGY_MEDOID_FRONTIER;

    /** attributes selected */
    BitSet attributes;

    /** m_Strategy: select the medoid */
    public static final int STRATEGY_MEDOID = 0;

    /** m_Strategy: select the medoid and frontier */
    public static final int STRATEGY_MEDOID_FRONTIER = 1;

    /** The m_Strategy used to select attributes */
    public static final Tag[] TAGS_STRATEGY = {
	new Tag(STRATEGY_MEDOID, "Select the medoid of each cluster"),
	new Tag(STRATEGY_MEDOID_FRONTIER, "Select the medoid and the frontier of each cluster"),};

    /**
     * Select the medoid of each cluster
     * @param kms partition obtained
     * @param data dataset partitioned
     * @param df distance function used
     * @throws Exception
     */
    protected void select_medoids(KMedoids kms, Instances data, SupervisedDistanceFunctionForSSF df) throws Exception {

	for (int k = 0; k < kms.numberOfClusters(); ++k) {
	    Instance medoid = kms.getClusterMedoids().instance(k);
	    attributes.set((int) medoid.value(0));
	}

    }

    /**
     * Select the frontier (feature farthest/less correlated with medoid and more correlated with the class)
     * of each cluster
     * @param kms partition obtained
     * @param data dataset partitioned
     * @param df distance function used
     * @throws Exception
     */
    protected void select_frontiers(KMedoidsSupervised kms, Instances data, SupervisedDistanceFunctionForSSF df) throws Exception {

	for (int k = 0; k < kms.numberOfClusters(); ++k) {
	    double greater = Double.NEGATIVE_INFINITY;
	    int greaterIndex = -1;
	    double dist;
	    Instance centroid = kms.getClusterMedoids().instance(k);
	    for (int i = 0; i < data.numInstances(); ++i) {
		if (kms.getAssignments()[i] != k) {
		    continue;
		}

		dist = df.distance(data.instance(i), centroid);

		//in supervised problems consider the class correlation
		dist = (dist + df.classCorrelation(data.instance(i))) / 2;


		if (dist > greater) {
		    greaterIndex = i;
		    greater = dist;
		}
	    }
	    //avoid singleton exception
	    if (greaterIndex != -1) {
		attributes.set(greaterIndex);
	    }
	}
    }

    /**
     * Initialize the evaluator
     * @param data dataset that will be used
     * @throws Exception
     */
    @Override
    public void buildEvaluator(Instances data) throws Exception {
	getCapabilities().testWithFail(data);
	/*
	if( data.classIndex() >= 0 && data.classIndex() != (data.numAttributes()-1) ){
	    throw new Exception("If set the class index must be the last attribute of the dataset");
	}
	 */
    }

    /**
     * Select features from the clusters induced
     * @param kms partition obtained
     * @param data dataset partitioned
     * @param df distance function used
     * @return a BitSet with the selected features
     * @throws Exception
     */
    public BitSet chooseFeatures(KMedoids kms, Instances data, DistanceFunctionForSSF df) throws Exception {
	if (!(df instanceof SupervisedDistanceFunctionForSSF)) {
	    throw new Exception("To run S3F a SupervisedDistanceFunction is needed");
	}
	if (!(kms instanceof KMedoidsSupervised)) {
	    throw new Exception("Use the supervised variant of KMedoids with S3F");
	}

	attributes = new BitSet(data.numAttributes());
	switch (m_Strategy) {
	    case STRATEGY_MEDOID:
		select_medoids((KMedoidsSupervised) kms, data, (SupervisedDistanceFunctionForSSF) df);
		break;
	    case STRATEGY_MEDOID_FRONTIER:
		select_medoids((KMedoidsSupervised) kms, data, (SupervisedDistanceFunctionForSSF) df);
		select_frontiers((KMedoidsSupervised) kms, data, (SupervisedDistanceFunctionForSSF) df);
		break;
	}

	return attributes;
    }

    /**
     * Set the m_Strategy to select features from the clusters
     * @param strat m_Strategy to use
     * @throws IllegalArgumentException
     */
    public void setStrategy(SelectedTag strat) throws IllegalArgumentException {
	if (strat.getTags() == TAGS_STRATEGY) {
	    m_Strategy = strat.getSelectedTag().getID();
	} else {
	    throw new IllegalArgumentException("Invalid strategy");
	}
    }

    /**
     * Return the m_Strategy to select features from the clusters
     * @return m_Strategy used
     */
    public SelectedTag getStrategy() {
	return new SelectedTag(m_Strategy, TAGS_STRATEGY);
    }

    /**
     * Tip text for this property
     * @return tip text for the GUI
     */
    public String strategyTipText() {
	return "The strategy to select features from the induced clusters";
    }

     /**
     * Returns an enumeration describing the available options.
     * @return an enumeration of all the available options.
     *
     **/
    public Enumeration<Option> listOptions() {
	Vector<Option> newVector = new Vector<Option>(1);

	newVector.addElement(new Option("\tSpecify the strategy to select features from induced clusters ("
		+ STRATEGY_MEDOID + " for medoids, and " + STRATEGY_MEDOID_FRONTIER + " for medoids and frontier).",
		"S", STRATEGY_MEDOID_FRONTIER, "-S <num>"));
	return newVector.elements();
    }

    /**
     * Parses a given list of options.
     *
    <!-- options-start -->
    * Valid options are: <p/>
    * 
    * <pre> -S &lt;num&gt;
    *  Specify the strategy to select features from induced clusters (0 for medoids, and 1 for medoids and frontier).</pre>
    * 
    <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    public void setOptions(String[] options) throws Exception {
	String optionString;
	resetOptions();

	optionString = Utils.getOption('S', options);
	if (optionString.length() != 0) {
	    setStrategy(new SelectedTag(Integer.parseInt(optionString), TAGS_STRATEGY));
	}
    }

    /**
     * Gets the current settings.
     * @return an array of strings suitable for passing to setOptions()
     */
    public String[] getOptions() {
	Vector<String> options = new Vector<String>();

	options.add("-S");
	options.add("" + m_Strategy);

	return options.toArray(new String[options.size()]);
    }

    /**
     * Reset options to default values
     */
    private void resetOptions() {
	m_Strategy = STRATEGY_MEDOID_FRONTIER;
    }

    /**
     * Returns default capabilities of the feature selection method.
     *
     * @return      the capabilities of this method
     */
    @Override
    public Capabilities getCapabilities() {
	Capabilities result = super.getCapabilities();
	result.disableAll();

	// attributes
	result.enable(Capability.NOMINAL_ATTRIBUTES);
	result.enable(Capability.NUMERIC_ATTRIBUTES);
	// class
	result.enable(Capability.NOMINAL_CLASS);
	

	return result;
    }

    /**
     * Main method for executing this evaluator.
     *
     * @param args the options, use "-h" to display options
     */
    public static void main(String[] args) {
	ASEvaluation.runEvaluator(new S3F(), args);
    }

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "S3F:\n\n"
		+ "Select features from an induced partition of the feature set considering " +
		"the correlation between features and class. "
		+ "See:\n\n" + getTechnicalInformation().toString();
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
	result.setValue(Field.AUTHOR, "Thiago F. Covões and Eduardo R. Hruschka");
	result.setValue(Field.YEAR, "2011");
	result.setValue(Field.TITLE, "Towards Improving Cluster-Based Feature Selection "
		+ "with a Simplified Silhouette Filter");
	result.setValue(Field.JOURNAL, "Information Sciences");
	result.setValue(Field.VOLUME, "181");
	result.setValue(Field.NUMBER, "18");
	result.setValue(Field.ISSN, "0020-0255");
	result.setValue(Field.PAGES, "3766-3782");

	return result;
    }

    /**
     * String representation
     */
    @Override
    public String toString(){
	return "";
    }
}

