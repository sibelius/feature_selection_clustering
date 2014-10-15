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
 *    InvertDataSet.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.preProcess;

import java.util.ArrayList;
import weka.core.Attribute;
import weka.core.Capabilities;
import weka.core.Instances;
import weka.core.Capabilities.Capability;
import weka.core.DenseInstance;
import weka.filters.SimpleBatchFilter;

/**
 * Filter used to transform features in instances
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public class InvertDataSet
	extends SimpleBatchFilter {

    /**
     * ID for serialization
     */
    private static final long serialVersionUID = 5650434109148256393L;
    /**
     * ID for instances (to identify which feature is each instance)
     */
    private int m_ID = 0;

    /**
     * Returns a string describing this filter
     *
     * @return a description of the filter suitable for
     * displaying in the explorer/experimenter gui
     */
    public String globalInfo() {
	return "Invert the dataset so the attributes becames the instances.";
    }

    /**
     * Returns the Capabilities of this filter.
     *
     * @return            the capabilities of this object
     * @see               Capabilities
     */
    @Override
    public Capabilities getCapabilities() {
	Capabilities result = super.getCapabilities();

	// attributes
	result.enableAllAttributes();
	result.enable(Capability.MISSING_VALUES);
	// class
	result.enableAllClasses();
	result.enable(Capability.NO_CLASS);
	return result;
    }

    /**
     * Define the format of the filtered dataset
     * @param input dataset to be filtered
     * @return format of the filtered dataset
     */
    protected Instances determineOutputFormat(Instances input) {
	ArrayList<Attribute> newAttributes = new ArrayList<Attribute>(input.numInstances() + 1);
	newAttributes.add(new Attribute("ID"));
	for (int i = 1; i <= input.numInstances(); ++i) {
	    newAttributes.add(new Attribute("Obj" + i));
	}
	int nAtributos = input.numAttributes();
	if (input.classIndex() >= 0) {
	    nAtributos--;
	}	
	Instances invertedData = new Instances(input.relationName(), newAttributes, nAtributos);
	setOutputFormat(invertedData);
	return invertedData;
    }

    /**
     * Filter a dataset
     * @param input dataset to be filtered
     * @return filtered dataset
     */
    protected Instances process(Instances input) {
	Instances invertedData = new Instances(determineOutputFormat(input), 0);
	int numAttributes = input.numAttributes();
	int numInstances = input.numInstances();

	for (int i = 0; i < numAttributes; ++i) {
	    if (input.classIndex() == i) {
		continue;
	    }
	    double[] vals = new double[numInstances + 1];
	    vals[0] = m_ID++;
	    for (int j = 0; j < numInstances; ++j) {
		vals[j + 1] = input.instance(j).value(i);
	    }
	    invertedData.add(new DenseInstance(1.0, vals));
	}
	return invertedData;
    }

    public static void main(String[] argv) {
	runFilter(new InvertDataSet(), argv);
    }
}
