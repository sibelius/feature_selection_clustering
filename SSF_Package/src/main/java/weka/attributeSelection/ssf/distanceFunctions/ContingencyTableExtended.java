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
 *    ContingencyTableExtended.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import weka.core.ContingencyTables;
import weka.core.Instance;
import weka.core.Instances;

/**
 * Helper class that generates contingency tables
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class ContingencyTableExtended extends ContingencyTables {

    /**
     * Generate a contingency table for use with
     * weka.core.ContingencyTable
     * @param data dataset to be considered
     * @param att1 index of one attribute
     * @param att2 index of one attribute
     * @return contingency table of att1 and att2
     */
    public static double[][] makeContingencyTable(Instances data, int att1, int att2) {

        int ni = data.attribute(att1).numValues();
        int nj = data.attribute(att2).numValues();

        double[][] counts = new double[ni][nj];

        for (int instIndex = 0; instIndex < data.numInstances(); ++instIndex) {
            Instance inst = data.instance(instIndex);
            counts[(int) inst.value(att1)][(int) inst.value(att2)]++;
        }

        return counts;
    }

    /**
     * Generate a contingency table with 3 features
     * @param data dataset to be considered
     * @param att1 index of one attribute
     * @param att2 index of one attribute
     * @param att3 index of one attribute
     * @return contingency table of the three features
     */
    public static double[][][] makeContingencyTable(Instances data, int att1, int att2, int att3) {

        int ni = data.attribute(att1).numValues();
        int nj = data.attribute(att2).numValues();
        int nk = data.attribute(att3).numValues();

        double[][][] counts = new double[ni][nj][nk];

        for (int instIndex = 0; instIndex < data.numInstances(); ++instIndex) {
            Instance inst = data.instance(instIndex);
            counts[(int) inst.value(att1)][(int) inst.value(att2)][(int) inst.value(att3)]++;
        }

        return counts;
    }
}
