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
 *    PreBuiltDistanceMatrix.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import weka.core.Instances;

/**
 * Class for distance functions that build the distance matrix
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public abstract class PreBuiltDistanceMatrix 
	extends DistanceFunctionForSSF{
	
	/**
	 * Instances that will be used 
	 */
	protected Instances m_Instances = null;
	
	/**
	 * Distance matrix
	 */
	protected double[][] m_DistMatrix = null;
	

	/**
	 * Sets the dataset being used
	 * @param insts dataset
	 */
	public void setInstances(Instances insts) {
		this.m_Instances = insts;
		try{
			if(m_DistMatrix == null){
				makeDistanceMatrix();
			}
		}catch(Exception e){
			e.printStackTrace();
		}
	
	}
	
	/**
	 * Clear the distance matrix
	 */
	public void resetMatrix(){
		m_DistMatrix =  null;
	}

	/**
	 * Generates the distance matrix
	 * @throws Exception
	 */
	protected abstract void makeDistanceMatrix() throws Exception;
}
