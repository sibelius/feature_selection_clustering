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
 *    DistanceFunctionForSSF.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import weka.core.Instance;
import weka.core.Instances;
import weka.core.OptionHandler;

/**
 * Class that represent distance functions that can be used in SSF.
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public abstract class DistanceFunctionForSSF implements OptionHandler {

    /**
     * Return the distance between two features.
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
