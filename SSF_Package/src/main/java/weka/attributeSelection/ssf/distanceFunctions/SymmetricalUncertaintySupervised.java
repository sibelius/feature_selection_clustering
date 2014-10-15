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
 *    SymmetricalUncertaintySupervised.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import java.io.Serializable;
import weka.core.ContingencyTables;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;
import weka.filters.Filter;

/**
<!-- globalinfo-start -->
* Symmetrical Uncertainty adapted to consider feature-class correlation. For more informatio see: Thiago F. Covões, Eduardo R. Hruschka (2011). Towards Improving Cluster-Based Feature Selection with a Simplified Silhouette Filter. Information Sciences. 181(18):3766-3782.
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
* <pre> -F &lt;discretizer spec&gt;
*  Specify the discretizer algorithm for continuous features</pre>
* 
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class SymmetricalUncertaintySupervised
	extends SymmetricalUncertainty
	implements DiscreteDistance, Serializable, SupervisedDistanceFunctionForSSF, TechnicalInformationHandler {

    private static final long serialVersionUID = -4602259179884514272L;
    private double[] m_corrAtributteClass;

    /**
     * Generate the distance matrix of a dataset
     * @throws Exception
     */
    @Override
    protected void makeDistanceMatrix() throws Exception {

	discFilter.setInputFormat(m_Instances);
	Instances cInstances = Filter.useFilter(m_Instances, discFilter);

	if (saveDiscreteData()) {
	    discreteData = cInstances;
	}

	int nAttributes = m_Instances.numAttributes();
	int classIndex = m_Instances.classIndex();

	m_DistMatrix = new double[nAttributes - 1][nAttributes - 1];

	m_corrAtributteClass = new double[nAttributes - 1];
	double[][] counts;
	int idxI, idxJ;
	for (int i = 0; i < nAttributes; ++i) {
	    if (i == classIndex) {
		continue;
	    }
	    //this makes the features that are after the class reduce their indexes
	    //by one as done by the invertDataset filter
	    idxI = i;
	    if (i > classIndex) {
		--idxI;
	    }
	    counts = ContingencyTableExtended.makeContingencyTable(cInstances, i, classIndex);
	    m_corrAtributteClass[idxI] = ContingencyTables.symmetricalUncertainty(counts);
	}

	for (int i = 0; i < nAttributes; ++i) {
	    if (i == classIndex) {
		continue;
	    }
	    //this makes the features that are after the class reduce their 
	    //indexes by one as done by the invertDataset filter
	    idxI = i;
	    if (idxI > classIndex) {
		--idxI;
	    }
	    //correlation between a feature and itself is maximum
	    m_DistMatrix[idxI][idxI] = 1;
	    for (int j = i + 1; j < nAttributes; ++j) {

		//if the class is in the middle of the features, skip it
		if(j == classIndex){
		    continue;
		}
		
		counts = ContingencyTableExtended.makeContingencyTable(cInstances, i, j);

		//same as above
		idxJ = j;
		if (idxJ > classIndex) {
		    --idxJ;
		}

		//SU stored in usual form to facilitate the distance computation (see below)
		m_DistMatrix[idxI][idxJ] =
			m_DistMatrix[idxJ][idxI] = ContingencyTables.symmetricalUncertainty(counts);

	    }

	}

    }

    /**
     * Return the distance between two features
     * @param first a feature inverted using the InvertDataset filter
     * @param second a feature inverted using the InvertDataset filter
     * @return the distance between the features
     */
    @Override
    public double distance(Instance first, Instance second) {
	int att1 = (int) first.value(0);
	int att2 = (int) second.value(0);
	double su = 1 - m_DistMatrix[att1][att2];
	su += Math.abs(m_corrAtributteClass[att1] - m_corrAtributteClass[att2]);
	return su / 2;
    }

    /**
     * Return the correlation between two features
     * @param first a feature inverted using the InvertDataset filter
     * @param second a feature inverted using the InvertDataset filter
     * @return the correlation between the features
     */
    public double correlation(Instance first, Instance second) {
	return m_DistMatrix[(int) first.value(0)][(int) second.value(0)];
    }

    /**
     * Return the correlation between a features and the class
     * @param inst a feature inverted using the InvertDataset filter
     * @return the distance between the feature inst and the class
     */
    public double classCorrelation(Instance inst) {
	return m_corrAtributteClass[(int) inst.value(0)];
    }

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    @Override
    public String globalInfo() {
	return "Symmetrical Uncertainty adapted to consider feature-class correlation."
		+ " For more informatio see: " + getTechnicalInformation().toString();
    }

    /**
     * Returns an instance of a TechnicalInformation object, containing
     * detailed information about the technical background of this class,
     * e.g., paper reference or book this class is based on.
     *
     * @return the technical information about this class
     */
    @Override
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

