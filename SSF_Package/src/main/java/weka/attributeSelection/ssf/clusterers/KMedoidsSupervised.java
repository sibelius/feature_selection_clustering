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
 *    KMedoidsSupervised.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.clusterers;

import weka.attributeSelection.ssf.distanceFunctions.SupervisedDistanceFunctionForSSF;
import weka.core.DenseInstance;
import weka.core.Instances;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;

/**
<!-- globalinfo-start -->
* Cluster data using the kMedoids algorithm with a modification in the medoid to consider feature-class correlation.See:<br/>
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
public class KMedoidsSupervised
	extends KMedoids
	implements TechnicalInformationHandler {

    /** for serialization */
    private static final long serialVersionUID = 4355242498335512227L;

     /**
     * Returns a string describing this clusterer
     * @return a description of the evaluator suitable for
     * displaying in the explorer/experimenter gui
     */
    @Override
    public String globalInfo() {
        return "Cluster data using the kMedoids algorithm with a modification in the medoid " +
		"to consider feature-class correlation." +
		"See:\n\n"+getTechnicalInformation().toString();
    }



    /**
     * Update the medoid of a cluster considering feature-feature correlation and
     * feature-class correlation
     * @param medoidIndex index of the medoid being updated
     * @param members instances on this cluster
     * @param updateClusterInfo flag to indicate update on the medoid vector
     * @return new medoids coordinates
     */
    @Override
    public double[] moveMedoid(int medoidIndex, Instances members, boolean updateClusterInfo) {
	double[] coords = new double[members.numAttributes()];
	double dist, minDist = Double.NEGATIVE_INFINITY;
	int medoid = -1;

	SupervisedDistanceFunctionForSSF su = (SupervisedDistanceFunctionForSSF) m_DistanceFunction;

	//singleton case
	if (members.numInstances() == 1) {
	    medoid = 0;
	} else {
	    for (int i = 0; i < members.numInstances(); ++i) {
		dist = 0;
		for (int j = 0; j < members.numInstances(); ++j) {
		    if (j == i) {
			continue;
		    }
		    dist += su.correlation(members.instance(i), members.instance(j));
		}
		dist = (dist / (members.numInstances() - 1)) + su.classCorrelation(members.instance(i));
		dist /= 2;

		if (dist > minDist) {
		    minDist = dist;
		    medoid = i;
		}

	    }
	}
	if (updateClusterInfo) {
	    Instances medoids = getClusterMedoids();
	    medoids.add(new DenseInstance(members.instance(medoid)));
	}
	return coords;
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

