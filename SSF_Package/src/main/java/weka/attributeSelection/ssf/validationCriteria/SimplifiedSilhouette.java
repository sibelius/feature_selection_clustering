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
 *    SimplifiedSilhouette.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.validationCriteria;

import java.io.Serializable;
import weka.attributeSelection.ssf.clusterers.KMedoids;
import weka.attributeSelection.ssf.distanceFunctions.DistanceFunctionForSSF;
import weka.core.Instances;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;

/**
<!-- globalinfo-start -->
* Simplified Silhouette relative criteria for partition quality assessment.For more information see:<br/>
* <br/>
* Eduardo R. Hruschka, Ricardo J.G.B. Campello, Leandro N. de Castro (2006). Evolving clusters in gene-expression data. Information Sciences. 176(13):1898-1927.
* <p/>
<!-- globalinfo-end -->
 *
<!-- technical-bibtex-start -->
* BibTeX:
* <pre>
* &#64;article{Hruschka2006,
*    author = {Eduardo R. Hruschka and Ricardo J.G.B. Campello and Leandro N. de Castro},
*    journal = {Information Sciences},
*    number = {13},
*    pages = {1898-1927},
*    title = {Evolving clusters in gene-expression data},
*    volume = {176},
*    year = {2006},
*    ISSN = {0020-0255}
* }
* </pre>
* <p/>
<!-- technical-bibtex-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class SimplifiedSilhouette
	extends ValidationCriteria
	implements Serializable, TechnicalInformationHandler {

    /** for serialization **/
    private static final long serialVersionUID = -8325841579780209386L;

    /**
     * Compute the quality of a partition
     * @param kms clusterer used
     * @param data data partitioned
     * @param df distance function used
     * @return average simplified silhouette of the partition
     */
    public double quality(KMedoids kms, Instances data, DistanceFunctionForSSF df) {
	double[] silh = new double[data.numInstances()];
	try {

	    double a, b;

	    //hold the distance between the instance and the centroids
	    double[] tmp = new double[kms.numberOfClusters()];
	    //cluster that each instance is assigned
	    int[] instanceAssignment = kms.getAssignments();

	    int[] clustersSize = kms.getClusterSizes();
	    Instances centroids = kms.getClusterMedoids();
	    for (int i = 0; i < data.numInstances(); ++i) {
		try {
		    //singleton case
		    if (clustersSize[instanceAssignment[i]] == 1) {
			silh[i] = 0;
			continue;
		    }
		} catch (Exception e) {
		    e.printStackTrace();

		}

		for (int k = 0; k < kms.numberOfClusters(); ++k) {
		    if (k == instanceAssignment[i]) {
			tmp[k] = Double.POSITIVE_INFINITY;
			continue;
		    }
		    tmp[k] = df.distance(data.instance(i), centroids.instance(k));
		}
		//intern distance
		a = df.distance(data.instance(i), centroids.instance(instanceAssignment[i]));

		//nearest neighbor distance (extern distance)
		b = tmp[Utils.minIndex(tmp)];
		silh[i] = (b - a) / (Math.max(a, b));
	    }


	} catch (Exception e) {
	    e.printStackTrace();
	    System.out.println("Problem calculating simplified silhouette."
		    + "\n" + e.getMessage());
	}

	return Utils.mean(silh);
    }

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "Simplified Silhouette relative criteria for partition quality assessment."
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
	result.setValue(Field.AUTHOR, "Eduardo R. Hruschka and Ricardo J.G.B. Campello and Leandro N. de Castro");
	result.setValue(Field.YEAR, "2006");
	result.setValue(Field.TITLE, "Evolving clusters in gene-expression data");
	result.setValue(Field.JOURNAL, "Information Sciences");
	result.setValue(Field.VOLUME, "176");
	result.setValue(Field.NUMBER, "13");
	result.setValue(Field.ISSN, "0020-0255");
	result.setValue(Field.PAGES, "1898-1927");

	return result;
    }
}

