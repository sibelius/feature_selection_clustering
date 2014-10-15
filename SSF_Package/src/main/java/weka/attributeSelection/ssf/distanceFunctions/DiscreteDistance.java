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
 *    DiscreteDistance.java
 *    Copyright (C) 2011 Thiago Covões
 *
 */
package weka.attributeSelection.ssf.distanceFunctions;

import weka.core.Instances;
import weka.core.OptionHandler;
import weka.filters.Filter;

/**
 * Interface for distance based on discrete features.
 * Define the method to set and get the discretized to use.
 *
 * @author Thiago Covões (tcovoes at icmc dot usp dot br)
 * @version 1.0
 */
public interface DiscreteDistance extends OptionHandler {

    /**
     * Sets the filter used to discretize the dataset
     * @param discretizer discretize filter
     */
    public void setDiscretizer(Filter discretizer);

    /**
     * Get the filter used to discretize the dataset
     * @return discretize filter
     */
    public Filter getDiscretizer();

    /**
     * Set the flag for saving the discrete dataset
     * @param save
     */
    public void saveDiscreteData(boolean save);

    /**
     * Get the flag which indicates if the discrete dataset is saved
     * @return true if discrete data is saved
     */
    public boolean saveDiscreteData();

    /**
     * Get the dataset discretized
     * @return the discretized dataset
     */
    public Instances getDiscreteData();
}
