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
 *    SupervisedDistanceFunctionForSSF.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import weka.core.Instance;
import weka.core.Instances;

/**
 * Interface to distance functions for SSF that consider feature-class correlation.
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public interface SupervisedDistanceFunctionForSSF {

    /**
     * Return the correlation of a feature with the class
     * @param inst a feature inverted using the InvertDataset filter
     * @return correlation between the feature and the class
     */
    public double classCorrelation(Instance inst);

    /**
     * Return the correlation between two features
     * @param obj1 a feature inverted using the InvertDataset filter
     * @param obj2 a feature inverted using the InvertDataset filter
     * @return correlation between features represented in obj1 and obj2
     */
    public double correlation(Instance obj1, Instance obj2);

    /**
     * Return the distance between two features
     * @param obj1 a feature inverted using the InvertDataset filter
     * @param obj2 a feature inverted using the InvertDataset filter
     * @return distance between features represented in obj1 and obj2
     */
    public abstract double distance(Instance obj1, Instance obj2);

    /**
     * Sets the dataset being used
     * @param data dataset
     */
    public abstract void setInstances(Instances data);
}
