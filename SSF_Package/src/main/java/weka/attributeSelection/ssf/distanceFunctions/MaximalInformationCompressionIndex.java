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
 *    MaximalInformationCompressionIndex.java
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
import weka.core.TechnicalInformation;
import weka.core.TechnicalInformation.Field;
import weka.core.TechnicalInformation.Type;
import weka.core.TechnicalInformationHandler;
import weka.core.Utils;

/**
<!-- globalinfo-start -->
* Maximal Information Compression Index. For more informatio see: Pabitra Mitra, C. A. Murthy, Sankar K. Pal (2002). Unsupervised Feature Selection Using Feature Similarity. IEEE Transactions on Pattern Analysis and Machine Intelligence. 24(3):301-312.
* <p/>
<!-- globalinfo-end -->
 *
<!-- technical-bibtex-start -->
* BibTeX:
* <pre>
* &#64;article{Mitra2002,
*    author = {Pabitra Mitra and C. A. Murthy and Sankar K. Pal},
*    journal = {IEEE Transactions on Pattern Analysis and Machine Intelligence},
*    number = {3},
*    pages = {301-312},
*    title = {Unsupervised Feature Selection Using Feature Similarity},
*    volume = {24},
*    year = {2002},
*    ISSN = {0162-8828}
* }
* </pre>
* <p/>
<!-- technical-bibtex-end -->
 *
<!-- options-start -->
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class MaximalInformationCompressionIndex
	extends PreBuiltDistanceMatrix
	implements DistanceForInvertedDataset, Serializable, TechnicalInformationHandler {

    /** for serialization */
    private static final long serialVersionUID = -7328528398438526450L;

    /**
     * Generate the distance matrix of a dataset
     *
     * @throws Exception
     */
    public void makeDistanceMatrix() throws Exception {

	int nInstances = m_Instances.numInstances();
	m_DistMatrix = new double[nInstances][nInstances];
	Instances cInstances = new Instances(m_Instances);
	cInstances.deleteAttributeAt(0);//remove id

	for (int i = 0; i < nInstances; ++i) {
	    double[] X = cInstances.instance(i).toDoubleArray();
	    double varX = Utils.variance(X);
	    //System.out.print("\n");
	    for (int j = i + 1; j < nInstances; ++j) {
		double[] Y = cInstances.instance(j).toDoubleArray();
		double varY = Utils.variance(Y);
		double corrCoef = Utils.correlation(X, Y, X.length);
		m_DistMatrix[i][j] =
			m_DistMatrix[j][i] = (varX + varY)
			- Math.sqrt(
			((varX + varY) * (varX + varY)) - 4 * varX * varY * (1 - (corrCoef * corrCoef)));
	    }
	}
    }

    /**
     * Get distance matrix
     * @return distance matrix
     */
    public double[][] getMatrix() {
	return m_DistMatrix;
    }

    /**
     * Return the distance between two features
     * @param first a feature inverted using the InvertDataset filter
     * @param second a feature inverted using the InvertDataset filter
     * @return distance between features represented in first and second
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
	return newVector.elements();
    }

    /**
     * Parses a given list of options.
     *
    <!-- options-start -->
    <!-- options-end -->
     *
     * @param options the list of options as an array of strings
     * @throws Exception if an option is not supported
     */
    public void setOptions(String[] options) throws Exception {
    }

    /**
     * Gets the current settings.
     * @return an array of strings suitable for passing to setOptions()
     */
    public String[] getOptions() {
	Vector<String> options = new Vector<String>();
	return options.toArray(new String[options.size()]);
    }

    /**
     * Returns a string describing this search method
     * @return a description of the search method suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "Maximal Information Compression Index. For more informatio see: " +
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

	result = new TechnicalInformation(Type.ARTICLE);
	result.setValue(Field.AUTHOR, "Pabitra Mitra and C. A. Murthy and Sankar K. Pal");
	result.setValue(Field.YEAR, "2002");
	result.setValue(Field.TITLE, "Unsupervised Feature Selection Using Feature Similarity");
	result.setValue(Field.JOURNAL, "IEEE Transactions on Pattern Analysis and Machine Intelligence");
	result.setValue(Field.VOLUME, "24");
	result.setValue(Field.NUMBER, "3");
	result.setValue(Field.ISSN, "0162-8828");
	result.setValue(Field.PAGES, "301-312");

	return result;
    }
}

