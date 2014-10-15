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
 *    KMedoids.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.clusterers;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Random;
import java.util.List;
import java.util.Vector;
import weka.attributeSelection.ssf.distanceFunctions.DistanceFunctionForSSF;
import weka.attributeSelection.ssf.distanceFunctions.PearsonCorrelation;

import weka.classifiers.rules.DecisionTableHashKey;
import weka.clusterers.NumberOfClustersRequestable;
import weka.clusterers.RandomizableClusterer;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.Utils;
import weka.core.WeightedInstancesHandler;
import weka.core.Capabilities.Capability;
import weka.core.DenseInstance;
import weka.core.OptionHandler;

/**
<!-- globalinfo-start -->
* Cluster data using the kMedoids algorithm.
* <p/>
<!-- globalinfo-end -->
 *
 *
 *
<!-- options-start -->
* Valid options are: <p/>
* 
* <pre> -N &lt;num&gt;
*  number of clusters.
*  (default 2).</pre>
* 
* <pre> -A &lt;classname and options&gt;
*  Distance function to use.
* </pre>
* 
* <pre> -I &lt;num&gt;
*  Maximum number of iterations.
* </pre>
* 
* <pre> -O
*  Preserve order of instances.
* </pre>
* 
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class KMedoids
        extends RandomizableClusterer
        implements NumberOfClustersRequestable, WeightedInstancesHandler, Serializable {

    /** for serialization */
    static final long serialVersionUID = -3235809600124455376L;

    /**
     * number of clusters to generate
     */
    private int m_NumClusters = 2;

    /**
     * holds the cluster medoids
     */
    private Instances m_ClusterMedoids;

    /**
     * The number of instances in each cluster
     */
    private int[] m_ClusterSizes;

    /**
     * Maximum number of iterations to be executed
     */
    private int m_MaxIterations = 10;

    /**
     * Keep track of the number of iterations completed before convergence
     */
    private int m_Iterations = 0;

    /**
     * Holds the squared errors for all clusters
     */
    private double[] m_squaredErrors;

    /** the distance function used. */
    protected DistanceFunctionForSSF m_DistanceFunction = new PearsonCorrelation();

    /**
     * Preserve order of instances
     */
    private boolean m_PreserveOrder = false;

    /**
     * Assignments obtained
     */
    protected int[] m_Assignments = null;

    /**
     * the default constructor
     */
    public KMedoids() {
        super();

        m_SeedDefault = 10;
        setSeed(m_SeedDefault);
    }

    /**
     * Returns a string describing this clusterer
     * @return a description of the evaluator suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
        return "Cluster data using the kMedoids algorithm.";
    }

    /**
     * Returns default capabilities of the clusterer.
     *
     * @return      the capabilities of this clusterer
     */
    @Override
    public Capabilities getCapabilities() {
        Capabilities result = super.getCapabilities();

        result.enable(Capability.NO_CLASS);
        // attributes
        result.enable(Capability.NOMINAL_ATTRIBUTES);
        result.enable(Capability.NUMERIC_ATTRIBUTES);
        result.enable(Capability.MISSING_VALUES);

        return result;
    }

    /**
     * Generates a clusterer. Has to initialize all fields of the clusterer
     * that are not being set via options.
     *
     * @param data set of instances serving as training data
     * @throws Exception if the clusterer has not been
     * generated successfully
     */
    public void buildClusterer(Instances data) throws Exception {

        // can clusterer handle the data?
        getCapabilities().testWithFail(data);

        m_Iterations = 0;

        Instances instances = new Instances(data);
  

        m_ClusterMedoids = new Instances(instances, m_NumClusters);
        int[] clusterAssignments = new int[instances.numInstances()];

        if (m_PreserveOrder) {
            m_Assignments = clusterAssignments;
        }

        m_DistanceFunction.setInstances(instances);


        initializeMedoids(instances);

        m_NumClusters = m_ClusterMedoids.numInstances();


        int i;
        boolean converged = false;
        int emptyClusterCount;
        Instances[] tempI = new Instances[m_NumClusters];
        m_squaredErrors = new double[m_NumClusters];

        while (!converged) {
            emptyClusterCount = 0;
            m_Iterations++;
            converged = true;
            for (i = 0; i < instances.numInstances(); i++) {
                Instance toCluster = instances.instance(i);
                int newC = clusterProcessedInstance(toCluster, true);
                if (newC != clusterAssignments[i]) {
                    converged = false;
                }
                clusterAssignments[i] = newC;
            }

            // update medoids
            m_ClusterMedoids = new Instances(instances, m_NumClusters);
            for (i = 0; i < m_NumClusters; i++) {
                tempI[i] = new Instances(instances, 0);
            }
            for (i = 0; i < instances.numInstances(); i++) {
                tempI[clusterAssignments[i]].add(instances.instance(i));
            }
            for (i = 0; i < m_NumClusters; i++) {
                if (tempI[i].numInstances() == 0) {
                    // empty cluster
                    emptyClusterCount++;
                } else {
                    moveMedoid(i, tempI[i], true);
                }
            }
            if (m_Iterations >= m_MaxIterations && emptyClusterCount == 0) {
                converged = true;
            }
            if (emptyClusterCount > 0) {
                m_NumClusters -= emptyClusterCount;
                if (converged) {
                    Instances[] t = new Instances[m_NumClusters];
                    int index = 0;
                    for (int k = 0; k < tempI.length; k++) {
                        if (tempI[k].numInstances() > 0) {
                            t[index++] = tempI[k];
                        }
                    }
                    tempI = t;
                } else {
                    tempI = new Instances[m_NumClusters];
                }
            }


            if (!converged) {
                m_squaredErrors = new double[m_NumClusters];

            }
        }


        m_ClusterSizes = new int[m_NumClusters];
        for (i = 0; i < m_NumClusters; i++) {
            m_ClusterSizes[i] = tempI[i].numInstances();
        }
    }

    /**
     * Move the medoid to it's new coordinates. Generate the medoid coordinates based
     * on it's  members (objects assigned to the cluster of the medoid) and the distance
     * function being used.
     * @param medoidIndex index of the medoid which the coordinates will be computed
     * @param members the objects that are assigned to the cluster of this medoid
     * @param updateClusterInfo if the method is supposed to update the m_Cluster arrays
     * @return the medoid coordinates
     */
    protected double[] moveMedoid(int medoidIndex, Instances members, boolean updateClusterInfo) {
        double[] coords = new double[members.numAttributes()];
        double dist, minDist = Double.MAX_VALUE;
        int medoid = -1;
        Instances medoids = getClusterMedoids();

        for (int i = 0; i < members.numInstances(); ++i) {
            dist = 0;
            for (int j = 0; j < members.numInstances(); ++j) {
                if (j == i) {
                    continue;
                }
                dist += m_DistanceFunction.distance(members.instance(i), members.instance(j));
            }

            if (dist < minDist) {
                minDist = dist;
                medoid = i;
            }

        }

        if (updateClusterInfo) {
            medoids.add(new DenseInstance(members.instance(medoid)));
        }
        return coords;
    }

    /**
     * clusters an instance that has been through the filters
     *
     * @param instance the instance to assign a cluster to
     * @param updateErrors if true, update the within clusters sum of errors
     * @return a cluster number
     */
    private int clusterProcessedInstance(Instance instance, boolean updateErrors) {
        double minDist = Integer.MAX_VALUE;
        int bestCluster = 0;
        for (int i = 0; i < m_NumClusters; i++) {
            double dist = m_DistanceFunction.distance(instance, m_ClusterMedoids.instance(i));
            if (dist < minDist) {
                minDist = dist;
                bestCluster = i;
            }
        }
        if (updateErrors) {
            m_squaredErrors[bestCluster] += minDist;
        }
        return bestCluster;
    }


    /**
     * Initialize the clusters according to the seed
     * @param instances dataset being clustered
     */
    private void initializeMedoids(Instances instances) throws Exception {
        Random RandomO = new Random(getSeed());
        int instIndex;
        HashMap<DecisionTableHashKey, Boolean> initC = new HashMap<DecisionTableHashKey, Boolean>();
        DecisionTableHashKey hk = null;

        List<Integer> indexesList = new ArrayList<Integer>(instances.numInstances());


        for (int j = 0; j < instances.numInstances(); j++) {
            indexesList.add(j);
        }

        Collections.shuffle(indexesList, RandomO);

        for (instIndex = instances.numInstances() - 1; instIndex >= 0; instIndex--) {

            hk = new DecisionTableHashKey(instances.instance(indexesList.get(instIndex)),
                    instances.numAttributes(), true);
            if (!initC.containsKey(hk)) {
                m_ClusterMedoids.add(instances.instance(indexesList.get(instIndex)));
                initC.put(hk, null);
            }
            if (m_ClusterMedoids.numInstances() == m_NumClusters) {
                break;
            }
        }
    }




    /**
     * Classifies a given instance.
     *
     * @param instance the instance to be assigned to a cluster
     * @return the number of the assigned cluster as an interger
     * if the class is enumerated, otherwise the predicted value
     * @throws Exception if instance could not be classified
     * successfully
     */
    @Override
    public int clusterInstance(Instance instance) throws Exception {
        return clusterProcessedInstance(instance, false);
    }

    /**
     * Returns the number of clusters.
     *
     * @return the number of clusters generated for a training dataset.
     * @throws Exception if number of clusters could not be returned
     * successfully
     */
    public int numberOfClusters() throws Exception {
        return m_NumClusters;
    }

    /**
     * Returns an enumeration describing the available options.
     *
     * @return an enumeration of all the available options.
     */
    @Override
    public Enumeration<Option> listOptions() {
        Vector<Option> result = new Vector<Option>();

        result.addElement(new Option(
                "\tnumber of clusters.\n"
                + "\t(default 2).",
                "N", 1, "-N <num>"));

        result.add(new Option(
                "\tDistance function to use.\n",
                "A", 1, "-A <classname and options>"));

        result.add(new Option(
                "\tMaximum number of iterations.\n",
                "I", 1, "-I <num>"));

        result.addElement(new Option(
                "\tPreserve order of instances.\n",
                "O", 0, "-O"));

        return result.elements();
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String numClustersTipText() {
        return "set number of clusters";
    }

    /**
     * set the number of clusters to generate
     *
     * @param n the number of clusters to generate
     * @throws Exception if number of clusters is negative
     */
    public void setNumClusters(int n) throws Exception {
        if (n <= 0) {
            throw new Exception("Number of clusters must be > 0");
        }
        m_NumClusters = n;
    }

    /**
     * gets the number of clusters to generate
     *
     * @return the number of clusters to generate
     */
    public int getNumClusters() {
        return m_NumClusters;
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String maxIterationsTipText() {
        return "set maximum number of iterations";
    }

    /**
     * set the maximum number of iterations to be executed
     *
     * @param n the maximum number of iterations
     * @throws Exception if maximum number of iteration is smaller than 1
     */
    public void setMaxIterations(int n) throws Exception {
        if (n <= 0) {
            throw new Exception("Maximum number of iterations must be > 0");
        }
        m_MaxIterations = n;
    }

    /**
     * gets the number of maximum iterations to be executed
     *
     * @return the number of clusters to generate
     */
    public int getMaxIterations() {
        return m_MaxIterations;
    }

    /**
     * Returns the tip text for this property.
     *
     * @return 		tip text for this property suitable for
     *         		displaying in the explorer/experimenter gui
     */
    public String distanceFunctionTipText() {
        return "The distance function to use for instances comparison";
    }

    /**
     * returns the distance function currently in use.
     *
     * @return the distance function
     */
    public DistanceFunctionForSSF getDistanceFunction() {
        return m_DistanceFunction;
    }

    /**
     * sets the distance function to use for instance comparison.
     *
     * @param df the new distance function to use
     * @throws Exception if instances cannot be processed
     */
    public void setDistanceFunction(DistanceFunctionForSSF df) throws Exception {
        m_DistanceFunction = df;
    }

    /**
     * Returns the tip text for this property
     * @return tip text for this property suitable for
     * displaying in the explorer/experimenter gui
     */
    public String preserveInstancesOrderTipText() {
        return "Preserve order of instances.";
    }

    /**
     * Sets whether order of instances must be preserved
     *
     * @param r true if missing values are to be
     * replaced
     */
    public void setPreserveInstancesOrder(boolean r) {
        m_PreserveOrder = r;
    }

    /**
     * Gets whether order of instances must be preserved
     *
     * @return true if missing values are to be
     * replaced
     */
    public boolean getPreserveInstancesOrder() {
        return m_PreserveOrder;
    }

    /**
     * Parses a given list of options. <p/>
     *
    <!-- options-start -->
    * Valid options are: <p/>
    * 
    * <pre> -N &lt;num&gt;
    *  number of clusters.
    *  (default 2).</pre>
    * 
    * <pre> -A &lt;classname and options&gt;
    *  Distance function to use.
    * </pre>
    * 
    * <pre> -I &lt;num&gt;
    *  Maximum number of iterations.
    * </pre>
    * 
    * <pre> -O
    *  Preserve order of instances.
    * </pre>
    * 
    <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    @Override
    public void setOptions(String[] options)
            throws Exception {

        String optionString = Utils.getOption('N', options);

        if (optionString.length() != 0) {
            setNumClusters(Integer.parseInt(optionString));
        }

        optionString = Utils.getOption("I", options);
        if (optionString.length() != 0) {
            setMaxIterations(Integer.parseInt(optionString));
        }

        String distFunctionClass = Utils.getOption('A', options);
        if (distFunctionClass.length() != 0) {
            String distFunctionClassSpec[] = Utils.splitOptions(distFunctionClass);
            if (distFunctionClassSpec.length == 0) {
                throw new Exception("Invalid DistanceFunction specification string.");
            }
            String className = distFunctionClassSpec[0];
            distFunctionClassSpec[0] = "";

            setDistanceFunction((DistanceFunctionForSSF) Utils.forName(DistanceFunctionForSSF.class,
                    className, distFunctionClassSpec));
        } else {
            setDistanceFunction(new PearsonCorrelation());
        }

        m_PreserveOrder = Utils.getFlag("O", options);

        super.setOptions(options);
    }

    /**
     * Gets the current settings
     *
     * @return an array of strings suitable for passing to setOptions()
     */
    @Override
    public String[] getOptions() {
        int i;
        Vector<String> result;
        String[] options;

        result = new Vector<String>();

        result.add("-N");
        result.add("" + getNumClusters());

        result.add("-A");
	if(m_DistanceFunction instanceof OptionHandler){
	    result.add((m_DistanceFunction.getClass().getName() + " "
		    + Utils.joinOptions(m_DistanceFunction.getOptions())).trim());
	}
        result.add("-I");
        result.add("" + getMaxIterations());

        if (m_PreserveOrder) {
            result.add("-O");
        }

        options = super.getOptions();
        for (i = 0; i < options.length; i++) {
            result.add(options[i]);
        }

        return result.toArray(new String[result.size()]);
    }

    /**
     * Reset options to default values
     */
    protected void resetOptions() {
	m_Assignments = null;
	m_ClusterMedoids = null;
	m_DistanceFunction = new PearsonCorrelation();
	m_Iterations = 0;
	m_MaxIterations = 10;
	m_NumClusters = 2;
	m_PreserveOrder = true;
    }


    private String pad(String source, String padChar,
            int length, boolean leftPad) {
        StringBuffer temp = new StringBuffer();

        if (leftPad) {
            for (int i = 0; i < length; i++) {
                temp.append(padChar);
            }
            temp.append(source);
        } else {
            temp.append(source);
            for (int i = 0; i < length; i++) {
                temp.append(padChar);
            }
        }
        return temp.toString();
    }

    /**
     * Gets the the cluster medoids
     *
     * @return		the cluster medoids
     */
    public Instances getClusterMedoids() {
        return m_ClusterMedoids;
    }

    /**
     * Gets the squared error for all clusters
     *
     * @return		the squared error
     */
    public double getSquaredError() {
        return Utils.sum(m_squaredErrors);
    }

    /**
     * Gets the number of instances in each cluster
     *
     * @return		The number of instances in each cluster
     */
    public int[] getClusterSizes() {
        return m_ClusterSizes;
    }

    /**
     * Gets the assignments for each instance
     * @return Array of indexes of the medoid assigned to each instance
     * @throws Exception if order of instances wasn't preserved or no assignments were made
     */
    public int[] getAssignments() throws Exception {
        if (!m_PreserveOrder) {
            throw new Exception("The assignments are only available when order of instances is preserved (-O)");
        }
        if (m_Assignments == null) {
            throw new Exception("No assignments made.");
        }
        return m_Assignments;
    }
}

