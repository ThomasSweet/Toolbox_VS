package de.sb.toolbox.util;

import de.sb.toolbox.Copyright;


/**
 * Instances of this interface are required to be immutable, and map objects to indices.
 */
@Copyright(year = 2009, holders = "Sascha Baumeister")
public interface Indexer {

	/**
	 * Returns the cardinality of the index calculation, i.e. the maximum index + 1.
	 * @return the cardinality, a strictly positive value
	 */
	int cardinality ();


	/**
	 * Returns the index calculated for the given object. The index returned must be guaranteed to be the same for equal values.
	 * @param object the object
	 * @return an index within range {@code [0, cardinality[}
	 * @throws NullPointerException if the given object cannot be indexed because it is {@code null}
	 * @throws IllegalArgumentException if the given object cannot be indexed for another reason
	 */
	int index (Object object) throws NullPointerException, IllegalArgumentException;
}