package de.sb.toolbox.xml;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import javax.xml.namespace.NamespaceContext;
import de.sb.toolbox.Copyright;


/**
 * Instances of this class model simple map based namespace contexts.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class MappedNamespaceContext implements NamespaceContext {

	private final Map<String,String> prefixes;


	/**
	 * Creates a new instance with no initial prefix registrations.
	 */
	public MappedNamespaceContext () {
		this.prefixes = new HashMap<>();
	}


	/**
	 * Creates a new instance with an initial prefix registration.
	 */
	public MappedNamespaceContext (final String prefix, final String namespaceURI) {
		this();
		if (prefix == null | namespaceURI == null) throw new NullPointerException();

		this.prefixes.put(prefix, namespaceURI);
	}


	/**
	 * Returns the prefixes. Note that the resulting map can be modified to change this context's prefix registrations.
	 * @return the prefix map
	 */
	public Map<String,String> getPrefixes () {
		return this.prefixes;
	}


	/**
	 * Returns the namespace URI registered for the given prefix.
	 * @param prefix the prefix
	 * @return the namespace URI
	 * @throws NoSuchElementException if the prefix is not registered within this context
	 */
	public String getNamespaceURI (final String prefix) {
		final String result = this.prefixes.get(prefix);
		if (result == null) throw new NoSuchElementException(prefix);
		return result;
	}


	/**
	 * Returns any prefix registered for the given namespace URI.
	 * @param namespaceURI the namespace URI
	 * @return a prefix
	 * @throws NoSuchElementException if the namespace URI is not registered within this context
	 */
	public String getPrefix (final String namespaceURI) {
		for (final Map.Entry<String,String> entry : this.prefixes.entrySet()) {
			if (namespaceURI.equals(entry.getValue())) return entry.getKey();
		}
		throw new NoSuchElementException(namespaceURI);
	}


	/**
	 * Returns all prefixes registered for the given namespace URI.
	 * @param namespaceURI the namespace URI
	 * @return the prefixes
	 */
	public Iterator<String> getPrefixes (final String namespaceURI) {
		final Set<String> prefixes = new HashSet<>();
		for (final Map.Entry<String,String> entry : this.prefixes.entrySet()) {
			if (namespaceURI.equals(entry.getValue())) prefixes.add(entry.getKey());
		}
		return prefixes.iterator();
	}
}