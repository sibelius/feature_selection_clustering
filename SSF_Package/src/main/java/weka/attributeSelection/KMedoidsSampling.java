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
 *    KMedoidsSampling.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection;

import weka.attributeSelection.ssf.distanceFunctions.DistanceFunctionForSSF;
import java.util.BitSet;
import java.util.Enumeration;
import java.util.Vector;

import weka.attributeSelection.ssf.preProcess.InvertDataSet;
import weka.attributeSelection.ssf.validationCriteria.ValidationCriteria;
import weka.attributeSelection.ssf.clusterers.KMedoids;
import weka.attributeSelection.ssf.clusterers.KMedoidsSupervised;
import weka.attributeSelection.ssf.distanceFunctions.DistanceForInvertedDataset;
import weka.attributeSelection.ssf.distanceFunctions.PreBuiltDistanceMatrix;
import weka.attributeSelection.ssf.distanceFunctions.SupervisedDistanceFunctionForSSF;
import weka.attributeSelection.ssf.distanceFunctions.SymmetricalUncertainty;
import weka.attributeSelection.ssf.validationCriteria.SimplifiedSilhouette;
import weka.core.Capabilities;
import weka.core.Capabilities.Capability;
import weka.core.CapabilitiesHandler;

import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;
import weka.filters.Filter;

/**
<!-- globalinfo-start -->
* kMedoids Sampling:<br/>
* <br/>
* Search for the best feature partition. To generate each partition the kMedoids algorithm is used varying the k parameter acoording to the defined options. For each value of k the kMedoids algorithm is executed a defined number of times.The partition quality is assessed by the validation criteria used.See:<br/>
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
* <pre> -P &lt;num&gt;
*  Specify the minimum number of attributes groups. (default=2)</pre>
* 
* <pre> -Q &lt;num&gt;
*  Specify the maximun number of attributes groups.</pre>
* 
* <pre> -N &lt;num&gt;
*  Number of repetitions for each value of K</pre>
* 
* <pre> -V &lt;location of the class&gt;
*  Validation criteria that will be used to assess partition quality. (default=Simplified Silhouette)</pre>
* 
* <pre> -D &lt;location of the class&gt;
*  Distance function that will be used. (default=Symmetrical Uncertainty)</pre>
* 
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class KMedoidsSampling
	extends ASSearch
	implements CapabilitiesHandler, OptionHandler, TechnicalInformationHandler {

    /**
     * for serialization
     */
    static final long serialVersionUID = 7841338639536821867L;
    /**
     * minimum number of clusters
     */
    protected int m_MinK = 2;
    /**
     * maximum number of clusters
     */
    protected int m_MaxK = 2;
    /**
     * validation criteria to be used to assess the partition quality
     */
    protected ValidationCriteria m_ValidationCriteria = new SimplifiedSilhouette();
    /**
     * number of repetitions to be executed for each value of k between m_minK and m_maxK
     */
    protected int m_NumberRepetitions = 20;
    /**
     * holds the merit of the best subset found
     */
    protected double bestMerit;
    /**
     * holds the k of the best subset found
     */
    protected double m_BestK;
    /**
     * best partition found during search
     */
    protected KMedoids bestPartition = null;
    /**
     * distance function to be used with kMedoids
     */
    protected DistanceFunctionForSSF m_DistanceFunction = new SymmetricalUncertainty();
    /**
     * Clusterer algorithm that will be used for partitioning the feature set
     */
    protected KMedoids m_Clusterer = new KMedoids();

    /**
     * Reference to the riginal dataset
     */
    protected Instances originalDataset = null;

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "kMedoids Sampling:\n\n"
		+ "Search for the best feature partition. To generate each partition the kMedoids algorithm is used "
		+ "varying the k parameter acoording to the defined options. "
		+ "For each value of k the kMedoids algorithm is executed a defined number of times."
		+ "The partition quality is assessed by the validation criteria used."
		+ "See:\n\n" + getTechnicalInformation().toString();
    }

    /**
     * KMedoids Sampling. Search the feature space through feature clustering.
     */
    public KMedoidsSampling() {
	resetOptions();
    }

    /**
     * Returns default capabilities of the feature selection method.
     *
     * @return      the capabilities of this method
     */
    public Capabilities getCapabilities() {
	Capabilities result = new Capabilities(this);
	result.disableAll();

	// attributes
	result.enable(Capability.NOMINAL_ATTRIBUTES);
	result.enable(Capability.NUMERIC_ATTRIBUTES);
	// class
	result.enable(Capability.NOMINAL_CLASS);
	result.enable(Capability.MISSING_CLASS_VALUES);
	//instances
	//actually this should be m_MinK but WEKA use this number to determine
	//the minimum number of instances with class value... and SSF does not need any class value
	result.setMinimumNumberInstances(0);
	return result;
    }

    /**
     * Returns an enumeration describing the available options.
     * @return an enumeration of all the available options.
     *
     **/
    public Enumeration<Option> listOptions() {
	Vector<Option> newVector = new Vector<Option>(5);

	newVector.addElement(new Option("\tSpecify the minimum number of attributes groups. (default=2)", "P", 1, "-P <num>"));
	newVector.addElement(new Option("\tSpecify the maximun number of attributes groups.", "Q", 1, "-Q <num>"));
	newVector.addElement(new Option("\tNumber of repetitions for each value of K", "N", 1, "-N <num>"));
	newVector.addElement(new Option("\tValidation criteria that will be used to assess partition quality."
		+ " (default=Simplified Silhouette)", "V", 1, "-V <location of the class>"));
	newVector.addElement(new Option("\tDistance function that will be used."
		+ " (default=Symmetrical Uncertainty)", "D", 1, "-D <location of the class>"));
	return newVector.elements();
    }

    /**
     * Parses a given list of options.
     *
    <!-- options-start -->
    * Valid options are: <p/>
    * 
    * <pre> -P &lt;num&gt;
    *  Specify the minimum number of attributes groups. (default=2)</pre>
    * 
    * <pre> -Q &lt;num&gt;
    *  Specify the maximun number of attributes groups.</pre>
    * 
    * <pre> -N &lt;num&gt;
    *  Number of repetitions for each value of K</pre>
    * 
    * <pre> -V &lt;location of the class&gt;
    *  Validation criteria that will be used to assess partition quality. (default=Simplified Silhouette)</pre>
    * 
    * <pre> -D &lt;location of the class&gt;
    *  Distance function that will be used. (default=Symmetrical Uncertainty)</pre>
    * 
    <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     *
     **/
    public void setOptions(String[] options)
	    throws Exception {
	String optionString;
	resetOptions();

	optionString = Utils.getOption('P', options);
	if (optionString.length() != 0) {
	    setMinK(Integer.parseInt(optionString));
	}

	optionString = Utils.getOption('Q', options);

	if (optionString.length() != 0) {
	    setMaxK(Integer.parseInt(optionString));
	}

	optionString = Utils.getOption('N', options);

	if (optionString.length() != 0) {
	    setNumberRepetitions(Integer.parseInt(optionString));
	}

	optionString = Utils.getOption('V', options);
	String [] validatorSpec = Utils.splitOptions(optionString);
        if (optionString.length() != 0) {
	    String validatorName = validatorSpec[0];
	    validatorSpec[0] = "";
	    setValidationCriteria((ValidationCriteria) Utils.forName(ValidationCriteria.class, validatorName, validatorSpec));
        }


	optionString = Utils.getOption('D', options);
	String [] distanceSpec = Utils.splitOptions(optionString);
	if (optionString.length() != 0) {
	    String distanceName = distanceSpec[0];
	    distanceSpec[0] = "";
	    setDistanceFunction((DistanceFunctionForSSF) Utils.forName(DistanceFunctionForSSF.class,
		    distanceName, distanceSpec));

	}
    }

    /**
     * Gets the current settings.
     * @return an array of strings suitable for passing to setOptions()
     */
    public String[] getOptions() {
	Vector<String> options = new Vector<String>();

	options.add("-P");
	options.add("" + m_MinK);
	options.add("-Q");
	options.add("" + m_MaxK);
	options.add("-N");
	options.add("" + m_NumberRepetitions);
	options.add("-V");
	if( m_ValidationCriteria instanceof OptionHandler){
	    options.add(m_ValidationCriteria.getClass().getName()+ " "
		    + Utils.joinOptions(((OptionHandler)m_ValidationCriteria).getOptions()));
	}else{
	    options.add("" + m_ValidationCriteria.getClass().getName());
	}
	options.add("-D");
	if( m_DistanceFunction instanceof OptionHandler){
	    options.add(m_DistanceFunction.getClass().getName()+ " "
		    + Utils.joinOptions(((OptionHandler)m_DistanceFunction).getOptions()));
	}else{
	    options.add("" + m_DistanceFunction.getClass().getName());
	}

	return options.toArray(new String[options.size()]);
    }

    /**
     * returns a description of the search as a String
     * @return a description of the search
     */
    @Override
    public String toString() {
	StringBuffer txt = new StringBuffer();
	txt.append("KMedoids Search.\nKmin: " + m_MinK + "\tKmax: " + m_MaxK
		+ "\tK: " + bestPartition.getNumClusters() + " \tBest Merit: " + bestMerit + "\n");
	txt.append("Clusterer used: " + getClusterer().getClass().getName() + "\n\n");
	try {
	    int[] assignments = bestPartition.getAssignments();
	    for (int c = 0; c < bestPartition.numberOfClusters(); ++c) {
		int attIdx =  (int) bestPartition.getClusterMedoids().instance(c).value(0);
		String attName = originalDataset.attribute( attIdx ).name();

		//in case the class was in the middle of the attributes
		//it is necessary to inc one from the idx to revert what is done by the invertdataset filter
		if(originalDataset.classIndex() >= 0 && attIdx >= originalDataset.classIndex()){
		    attName = originalDataset.attribute( attIdx+1 ).name();
		}
		txt.append("Cluster " + (c + 1)
			+ ":\tMedoid: " + attName + "\n\n");
		for (int i = 0; i < assignments.length; ++i) {
		    if (assignments[i] == c) {
			attName = originalDataset.attribute(i).name();
			//in case the class was in the middle of the attributes as above
			if(originalDataset.classIndex() >= 0 && i >= originalDataset.classIndex()){
			   attName = originalDataset.attribute( i+1 ).name();
			}
			txt.append(attName + " , ");
		    }
		}
		txt.append("\n\n");
	    }

	} catch (Exception e) {
	    e.printStackTrace();
	}
	return txt.toString();

    }

    /**
     * Searches the attribute subset space by kMedoids Sampling search
     *
     * @param ASEval the attribute picker, must be a SSF
     * @param data the training instances.
     * @return an array (not necessarily ordered) of selected attribute indexes
     * @throws Exception if the search can't be completed
     */
    public int[] search(ASEvaluation ASEval, Instances data)
	    throws Exception {

	getCapabilities().testWithFail(data);

	if (!(ASEval instanceof SSFInterface)) {
	    throw new Exception(ASEval.getClass().getName()
		    + " is not a "
		    + "SSF!");
	}

	if (!(m_DistanceFunction instanceof SupervisedDistanceFunctionForSSF)
		&& (ASEval instanceof S3F)) {
	    throw new Exception("To run S3F a SupervisedDistanceFunction must be used");
	}

	/*
	if( data.classIndex() >= 0 && data.classIndex() != (data.numAttributes()-1) ){
	    throw new Exception("If set the class index must be the last attribute of the dataset");
	}
	 */
	if (m_DistanceFunction instanceof SupervisedDistanceFunctionForSSF) {
	    setClusterer(new KMedoidsSupervised());
	}

	if (m_DistanceFunction instanceof PreBuiltDistanceMatrix) {
	    ((PreBuiltDistanceMatrix) m_DistanceFunction).resetMatrix();
	}

	if (!(m_DistanceFunction instanceof DistanceForInvertedDataset)) {
	    m_DistanceFunction.setInstances(data);
	}

	int classIndex = data.classIndex();
	originalDataset = data;
	InvertDataSet invertData = new InvertDataSet();
	invertData.setInputFormat(data);
	Instances invertedDataSet = Filter.useFilter(data, invertData);

	if (m_DistanceFunction instanceof DistanceForInvertedDataset) {
	    m_DistanceFunction.setInstances(invertedDataSet);
	}

	SSFInterface ASEvaluator = (SSFInterface) ASEval;

	KMedoids kms = null;
	double obtained = 0;
	bestMerit = Double.NEGATIVE_INFINITY;


	for (int k = m_MinK; k <= m_MaxK; ++k) {

	    for (int rpt = 0; rpt < m_NumberRepetitions; ++rpt) {

		kms = getClusterer();
		KMedoids kmsExec = (KMedoids) KMedoids.makeCopy(kms);
		kmsExec.setDistanceFunction(m_DistanceFunction);
		kmsExec.setPreserveInstancesOrder(true);
		kmsExec.setNumClusters(k);
		kmsExec.setSeed(k * (rpt + 1));
		kmsExec.setMaxIterations(50);
		kmsExec.buildClusterer(invertedDataSet);
		obtained = m_ValidationCriteria.quality(kmsExec, invertedDataSet, m_DistanceFunction);

		if (Utils.gr(obtained, bestMerit)) {

		    bestMerit = obtained;
		    bestPartition = kmsExec;
		    m_BestK = kmsExec.getNumClusters();
		}
	    }

	}

	BitSet attributes = ASEvaluator.chooseFeatures(bestPartition, invertedDataSet, m_DistanceFunction);

	//release memory used by distance function to gc
	if (m_DistanceFunction instanceof PreBuiltDistanceMatrix) {
	    ((PreBuiltDistanceMatrix) m_DistanceFunction).resetMatrix();
	}

	int []attrs_selected = attributeList(attributes, invertedDataSet.numInstances());
	//fix for the case that the class attribute is not the last
	//as the class is removed by the invertDataset filter the attributes after it
	//have their indices decremented by one, the following code just fix this
	for(int i=0;i<attrs_selected.length;++i){
	    if(attrs_selected[i] >= classIndex){
		++attrs_selected[i];
	    }
	}

	return attrs_selected;
    }

    /**
     * Reset options to default values
     */
    protected void resetOptions() {
	bestMerit = 0;
	originalDataset = null;
	m_DistanceFunction = new SymmetricalUncertainty();
	m_ValidationCriteria = new SimplifiedSilhouette();
	m_MinK = 2;
	m_MaxK = 2;
	m_NumberRepetitions = 2;
    }

    /**
     * converts a BitSet into a list of attribute indexes
     * @param group the BitSet to convert
     * @param numAttribs
     * @return an array of attribute indexes
     **/
    protected int[] attributeList(BitSet group, int numAttribs) {
	int count = 0;

	// count how many were selected
	for (int i = 0; i < numAttribs; i++) {
	    if (group.get(i)) {
		count++;
	    }
	}

	int[] list = new int[count];
	count = 0;

	for (int i = 0; i < numAttribs; i++) {
	    if (group.get(i)) {
		list[count++] = i;
	    }
	}

	return list;
    }

    /**
     * Set the clustering algorithm to use
     * @param clusterer clustering algorithm
     */
    public void setClusterer(KMedoids clusterer) {
	this.m_Clusterer = clusterer;
    }

    /**
     * Gets the m_Clusterer algorithms being used
     * @return the m_Clusterer
     */
    public KMedoids getClusterer() {
	return this.m_Clusterer;
    }

    public String clustererTipText() {
	return "Clusterer that will be used to cluster features";
    }

    /**
     * Gets the distance function being used
     * @return the m_DistanceFunction
     */
    public DistanceFunctionForSSF getDistanceFunction() {
	return m_DistanceFunction;
    }

    /**
     * Set the distance function that will be used
     * @param function the distance function to set
     */
    public void setDistanceFunction(DistanceFunctionForSSF function) {
	m_DistanceFunction = function;
    }

    /**
     * Returns the tip text for this property.
     * @return tip text for this property suitable for displaying in the explorer/experimenter gui
     */
    public String distanceFunctionTipText() {
	return "The distance function that will be used to measure dissimilarity or "
		+ "the lack of correlation";
    }

    /**
     * Get the number of different initial partitions for each number of clusters
     * @return the m_NumberRepetitions
     */
    public int getNumberRepetitions() {
	return m_NumberRepetitions;
    }

    /**
     * Set the number of different initial partition for each number of clusters
     * @param repetitions the number of partitions to set
     */
    public void setNumberRepetitions(int repetitions) {
	m_NumberRepetitions = repetitions;
    }

    /**
     * Returns the tip text for this property.
     * @return tip text for this property suitable for displaying in the explorer/experimenter gui
     */
    public String numberRepetitionsTipText() {
	return "Number of different initial partitions for kMedoids for each number of clusters";
    }

    /**
     * Get the relative validation criteria being used
     * @return the m_ValidationCriteria
     */
    public ValidationCriteria getValidationCriteria() {
	return m_ValidationCriteria;
    }

    /**
     * Set the relative validation criteria to use
     * @param criteria the validation criteria to set
     */
    public void setValidationCriteria(ValidationCriteria criteria) {
	m_ValidationCriteria = criteria;
    }

    /**
     * Returns the tip text for this property.
     * @return tip text for this property suitable for displaying in the explorer/experimenter gui
     */
    public String validationCriteriaTipText() {
	return "Relative validation criteria that will be used to select the best partition";
    }

    /**
     * Get the minimum number of clusters
     * @return the m_MinK
     */
    public int getMinK() {
	return m_MinK;
    }

    /**
     * Set the minimum number of clusters
     * @param m the minimum number of clusters
     */
    public void setMinK(int m) {
	this.m_MinK = m;
    }

    /**
     * Returns the tip text for this property.
     * @return tip text for this property suitable for displaying in the explorer/experimenter gui
     */
    public String minKTipText() {
	return "Minimum number of clusters";
    }

    /**
     * Set the maximum number of clusters
     * @param m the maximum number of clusters
     */
    public void setMaxK(int m) {
	this.m_MaxK = m;
    }

    /**
     * Get the maximum number of clusters
     * @return the m_MaxK
     */
    public int getMaxK() {
	return this.m_MaxK;
    }

    /**
     * Returns the tip text for this property.
     * @return tip text for this property suitable for displaying in the explorer/experimenter gui
     */
    public String maxKTipText() {
	return "Maximun number of clusters";
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
}

