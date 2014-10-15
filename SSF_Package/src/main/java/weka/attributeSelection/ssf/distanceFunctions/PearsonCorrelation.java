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
 *    PearsonCorrelation.java
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
import weka.core.Utils;

/**
<!-- globalinfo-start -->
* The well known Pearson Correlation
* <p/>
<!-- globalinfo-end -->
 *
 *
<!-- options-start -->
<!-- options-end -->
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class PearsonCorrelation
	extends PreBuiltDistanceMatrix
	implements DistanceForInvertedDataset, Serializable {

    private static final long serialVersionUID = 3707730296838115813L;

    /**
     * Generate the distance matrix of a dataset
     * @throws Exception
     */
    public void makeDistanceMatrix() throws Exception {
	int nInstances = m_Instances.numInstances();
	m_DistMatrix = new double[nInstances][nInstances];
	Instances cInstances = new Instances(m_Instances);
	cInstances.deleteAttributeAt(0);//remove id

	for (int i = 0; i < nInstances; ++i) {
	    double[] X = cInstances.instance(i).toDoubleArray();
	    for (int j = i + 1; j < nInstances; ++j) {
		double[] Y = cInstances.instance(j).toDoubleArray();
		m_DistMatrix[i][j] =
			m_DistMatrix[j][i] = 1 - Math.abs(Utils.correlation(X, Y, X.length));
	    }
	}
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
     *
     **/
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
	return "The well known Pearson Correlation";
    }
}

