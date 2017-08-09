package de.sb.toolbox.net;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Executable;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import javax.xml.bind.annotation.adapters.XmlAdapter;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.util.StringEntry;


/**
 * Some Java types do not map naturally to a XML representation, for example {@Code HashMap} or other non JavaBean classes.
 * Conversely, a XML repsentation may map to a Java type but an application may choose to accesss the XML representation using
 * another Java type. This facade provides adapter classes adapting a bound type to a value type or vice versa. Their methods
 * are invoked by the JAXB binding framework during marshaling and unmarshalling
 */
@Copyright(year = 2012, holders = "Sascha Baumeister")
public final class XmlAdapters {

	/**
	 * Prevents external instantiation.
	 */
	private XmlAdapters () {}



	/**
	 * JAB-B adapter that marshales Path instances into byte[], and vice versa. This implies that the file content is marshaled,
	 * not the file path itself.
	 */
	static public class PathToBytes extends XmlAdapter<byte[],Path> {

		/**
		 * {@inheritDoc}
		 * @throws IOException if the given path doesn't represent a readable file, or if there is an I/O related problem
		 */
		@Override
		public byte[] marshal (final Path path) throws IOException {
			return path == null ? null : Files.readAllBytes(path);
		}


		/**
		 * {@inheritDoc}
		 * @throws IOException if there is an I/O related problem
		 */
		@Override
		public Path unmarshal (final byte[] content) throws IOException {
			if (content == null) return null;
			final Path path = Files.createTempFile("temp-", ".bin");
			Files.write(path, content);
			return path;
		}
	}



	/**
	 * JAB-B adapter that marshales Path instances into String, and vice versa. This implies that the file path is marshaled,
	 * not the file content.
	 */
	public class PathToString extends XmlAdapter<String,Path> {

		/**
		 * {@inheritDoc}
		 */
		@Override
		public String marshal (final Path path) {
			return path == null ? null : path.toString();
		}


		/**
		 * {@inheritDoc}
		 * @throws InvalidPathException if the given text cannot be converted to a {@code Path}
		 */
		@Override
		public Path unmarshal (final String text) throws InvalidPathException {
			return text == null ? null : Paths.get(text);
		}
	}



	/**
	 * JAB-B adapter that marshales Object instances into String using {@link Object#toString()}, and vice versa using a static
	 * {@code valueOf(String)} method, or a constructor featuring a text argument. This works for all primitive type wrappers
	 * (Character has special support), String and related types, enums, ...
	 */
	static public class ObjectToString extends XmlAdapter<String,Object> {

		/**
		 * {@inheritDoc}
		 */
		@Override
		public String marshal (final Object value) throws IllegalArgumentException {
			if (value == null) return null;
			if (value instanceof Character) return value.toString();

			stringFactory(value.getClass());// check unmarshal support
			return value.getClass().getName() + "@" + value.toString();
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public Object unmarshal (final String value) throws IllegalArgumentException {
			if (value == null) return null;
			if (value.length() == 1) return value.charAt(0);

			final int delimiterOffset = value.indexOf('@');
			final String typeName = value.substring(0, delimiterOffset);
			final String textValue = value.substring(delimiterOffset + 1);

			try {
				final Class<?> type = Class.forName(typeName, true, Thread.currentThread().getContextClassLoader());

				final Executable factory = stringFactory(type);
				if (factory instanceof Method) return ((Method) factory).invoke(null, textValue);
				if (factory instanceof Constructor) return ((Constructor<?>) factory).newInstance(textValue);

				throw new AssertionError();
			} catch (final ClassNotFoundException | InstantiationException | IllegalAccessException exception) {
				throw new IllegalArgumentException(exception);
			} catch (final InvocationTargetException exception) {
				throw new IllegalArgumentException(exception.getCause());
			}
		}


		/**
		 * Returns the factory executable used for unmarshaling.
		 * @param type the type
		 * @return the factory executable
		 * @throws NullPointerException if the given argument is {@code null}
		 * @throws IllegalArgumentException if the given argument is illegal
		 */
		static private Executable stringFactory (final Class<?> type) throws NullPointerException, IllegalArgumentException {
			if (type == String.class) {
				try {
					return type.getDeclaredMethod("valueOf", Object.class);
				} catch (final NoSuchMethodException exception) {
					throw new AssertionError(exception);
				}
			}

			for (Class<?> valueType = type; valueType != null; valueType = valueType.getSuperclass()) {
				try {
					final Method method = valueType.getDeclaredMethod("valueOf", String.class);
					if (method.getReturnType().isPrimitive()) throw new IllegalArgumentException();
					method.setAccessible(true);
					return method;
				} catch (final NoSuchMethodException exception) {
					// continue
				}
			}

			try {
				final Constructor<?> constructor = type.getDeclaredConstructor(String.class);
				constructor.setAccessible(true);
				return constructor;
			} catch (final NoSuchMethodException exception) {
				throw new IllegalArgumentException();
			}
		}
	}



	/**
	 * JAB-B adapter that marshales string maps to arrays using both {@link StringEntry} and {@link ObjectToString}, and vice
	 * versa. This works for all values that are primitive type wrappers (Character has special support), String and related
	 * types, enums, ...
	 */
	static public class MapToArray extends XmlAdapter<StringEntry[],Map<String,?>> {
		private final ObjectToString delegate = new ObjectToString();


		@Override
		public StringEntry[] marshal (final Map<String,?> map) {
			if (map == null) return null;

			@SuppressWarnings({ "unchecked", "rawtypes" })
			final Iterator<Entry<String,?>> iterator = (Iterator) map.entrySet().iterator();
			final StringEntry[] array = new StringEntry[map.size()];
			for (int index = 0; index < array.length; ++index) {
				final Entry<String,?> entry = iterator.next();
				final String key = entry.getKey();
				final String value = this.delegate.marshal(entry.getValue());
				array[index] = new StringEntry(key, value);
			}

			return array;
		}


		@Override
		public Map<String,?> unmarshal (final StringEntry[] array) {
			if (array == null) return null;

			final Map<String,Object> map = new HashMap<>();
			for (final StringEntry entry : array) {
				final String key = entry.getKey();
				final Object value = this.delegate.unmarshal(entry.getValue());
				map.put(key, value);
			}
			return map;
		}
	}
}