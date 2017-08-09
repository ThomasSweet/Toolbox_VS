package de.sb.toolbox.xml;

import java.io.StringWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.regex.Pattern;
import javax.xml.namespace.QName;
import de.sb.toolbox.Copyright;


/**
 * This facade provides methods to convert Java classes into/from equivalent QNames, and Java packages into/from equivalent
 * URIs. The mapping interprets a package name as a canonical host name in reversed component order, and vice versa. Note that
 * this class is declared final because it is a facade, and therefore not supposed to be extended.
 */
@Copyright(year = 2010, holders = "Sascha Baumeister")
public final class Namespaces {
	static private final Pattern DOT_PATTERN = Pattern.compile("\\.");


	/**
	 * Prevents instantiation.
	 */
	private Namespaces () {}


	/**
	 * Splits the given test into components divided by dots, reverses the component's order, reassembles the reversed
	 * components into a dot-divided string, and returns the result.
	 * @param text the dot-divided text to be reversed
	 * @return the reversed dot-divided text
	 * @throws NullPointerException if the given text is {@code null}
	 */
	static private String reverseComponents (final String text) {
		final String[] components = DOT_PATTERN.split(text);
		if (components.length < 2) return text;

		final StringWriter charSink = new StringWriter();
		charSink.write(components[components.length - 1]);
		for (int index = components.length - 2; index >= 0; --index) {
			charSink.write('.');
			charSink.write(components[index]);
		}

		return charSink.toString();
	}


	/**
	 * Returns a qualified name for the given class
	 * @param clazz the class
	 * @return the qualified name
	 * @throws NullPointerException if the given class is {@code null}
	 */
	static public QName toQualifiedName (final Class<?> clazz) {
		final URI packageNamespace = Namespaces.toURI(clazz.getPackage());
		return new QName(packageNamespace.toASCIIString(), clazz.getSimpleName());
	}


	/**
	 * Returns a class for the given qualified name
	 * @param qualifiedName the qualified name
	 * @return the class
	 * @throws NullPointerException if the given name is {@code null}
	 * @throws IllegalArgumentException if the given name doesn't represent a valid class
	 */
	static public Class<?> toClass (final QName qualifiedName) {
		if (qualifiedName == null) throw new NullPointerException();

		try {
			final URI packageNamespace = new URI(qualifiedName.getNamespaceURI());
			final String className = reverseComponents(packageNamespace.getHost()) + "." + qualifiedName.getLocalPart();
			return Class.forName(className, true, Thread.currentThread().getContextClassLoader());
		} catch (final Exception exception) {
			throw new IllegalArgumentException(exception);
		}
	}


	/**
	 * If the given package is {@code null} or empty, {@code http:///} is returned to represent Java's virtual default package.
	 * Otherwise returns a URI containing the given scheme, the component-reversed package name as host name, and {@code "/"} as
	 * path.
	 * @param pkg the package
	 * @return the package's equivalent URI, or {@code [scheme:]///} representing the default package
	 * @throws NullPointerException if the given package is {@code null}
	 * @throws IllegalArgumentException if the given package cannot be converted into a valid URI
	 */
	static public URI toURI (final Package pkg) {
		try {
			return new URI("http", pkg == null ? null : reverseComponents(pkg.getName()), "/", null);
		} catch (final URISyntaxException exception) {
			throw new IllegalArgumentException(pkg.toString());
		}
	}


	/**
	 * If the given URI's is {@code http:/}, {@code null} is returned to represent Java's virtual default package. Otherwise,
	 * this method returns a Java package that matches the component-reversed host of the given URI.
	 * @param uri the URI
	 * @return the URI's equivalent package, or {@code null} representing the default package
	 * @throws NullPointerException if the given URI is {@code null}
	 * @throws IllegalArgumentException if the given URI's scheme is not {@code "http"}, it's port or path are not undefined, or
	 *         it's host cannot be converted into a valid package accessible by the current thread's class loader
	 * @see QName
	 */
	static public Package toPackage (final URI uri) {
		if (!uri.getScheme().equals("http")) throw new IllegalArgumentException();
		if (uri.getPort() != -1) throw new IllegalArgumentException();
		if (uri.getPath() != null & !uri.getPath().equals("/")) throw new IllegalArgumentException();
		if (uri.getHost() == null || uri.getHost().isEmpty()) return null;

		final Package pkg = Package.getPackage(reverseComponents(uri.getHost()));
		if (pkg == null) throw new IllegalArgumentException(uri.toString());
		return pkg;
	}
}