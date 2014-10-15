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
 *    SSFInterface.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection;


import java.util.BitSet;

import weka.attributeSelection.ssf.clusterers.KMedoids;
import weka.attributeSelection.ssf.distanceFunctions.DistanceFunctionForSSF;
import weka.core.Instances;

/**
 * Interface for implementations of SSF feature selectors
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public interface SSFInterface {

    /**
     * Select features from the clusters induced
     * @param kms partition obtained
     * @param data dataset partitioned
     * @param df distance function used
     * @return a BitSet with the selected features
     * @throws Exception
     */
    public BitSet chooseFeatures(KMedoids kms, Instances data, DistanceFunctionForSSF df) throws Exception;
       
       
}
