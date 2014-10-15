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
 *    ValidationCriteria.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.validationCriteria;

import weka.attributeSelection.ssf.clusterers.KMedoids;
import weka.attributeSelection.ssf.distanceFunctions.DistanceFunctionForSSF;
import weka.core.Instances;

/**
 * Abstract class for relative validation criteria to be used for choosing the best partition in SSF
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public abstract class ValidationCriteria {

    /**
     * Check the quality of the obtained partition
     * @param kms kMedoids instance used in the clustering process
     * @param data Clustered data
     * @param df Distance function that will be used
     * @return Partition quality
     */
    abstract public double quality(KMedoids kms, Instances data, DistanceFunctionForSSF df);
}
