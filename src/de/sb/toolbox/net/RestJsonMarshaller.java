package de.sb.toolbox.net;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.StringWriter;
import java.io.Writer;
import java.lang.annotation.Annotation;
import java.lang.reflect.Array;
import java.lang.reflect.Type;
import java.net.ProtocolException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import javax.ws.rs.Consumes;
import javax.ws.rs.Produces;
import javax.ws.rs.WebApplicationException;
import javax.ws.rs.core.Feature;
import javax.ws.rs.core.FeatureContext;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.MultivaluedMap;
import javax.ws.rs.ext.MessageBodyReader;
import javax.ws.rs.ext.MessageBodyWriter;
import javax.ws.rs.ext.Provider;
import de.sb.toolbox.Copyright;


/**
 * This feature provides support for generic JSON marshalling in JAX-RS based REST servers.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class RestJsonMarshaller implements Feature {

	/**
	 * {@inheritDoc}
	 */
	public boolean configure (final FeatureContext context) {
		context.register(GenericJsonEntityWriter.class);
		context.register(GenericJsonEntityReader.class);
		context.register(GenericJsonListReader.class);
		return true;
	}



	/**
	 * Instances of this class marshal generic objects into JSON representations. A given object must be a primitive value, a
	 * string, or either an array, a collection or a string-map recursively containing the aforementioned types. Note that this
	 * class is designed to marshal raw JDBC table data into JSON, not JAX-B annotated entities (use MOXY instead). Also note
	 * that this provider does not actively check for cyclic objects, so any recursive entity structure will cause a stack
	 * overflow during marshalling.
	 */
	@Provider
	@Produces(MediaType.APPLICATION_JSON)
	static private class GenericJsonEntityWriter implements MessageBodyWriter<Object> {

		/**
		 * {@inheritDoc}
		 */
		@Override
		public boolean isWriteable (final Class<?> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType) {
			return MediaType.APPLICATION_JSON_TYPE.isCompatible(mediaType);
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public long getSize (final Object entity, final Class<?> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType) {
			return -1;
		}


		/**
		 * Writes the given entity to the given writer in JSON representation. Note that data is written into a temporary stream
		 * first, and only copied into the given entity stream once marshaling completes successfully. This is required to
		 * retain the ability to set response codes during exception handling, i.e. before the HTTP header content is committed.
		 * Note that UTF-8 encoding is used for all data.
		 * @throws NullPointerException if the given entity, media type or entity stream is {@code null}
		 * @throws IllegalArgumentException if the given media type is not JSON
		 * @throws IOException if there is an I/O related problem
		 */
		@Override
		public void writeTo (final Object entity, final Class<?> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType, final MultivaluedMap<String,Object> httpHeaders, final OutputStream entityStream) throws NullPointerException, IOException {
			assert this.isWriteable(type, genericType, annotations, mediaType);
			try (Writer entityWriter = new StringWriter()) {
				RestJsonMarshaller.marshal(entity, entityWriter);

				final Writer streamWriter = new OutputStreamWriter(entityStream, StandardCharsets.UTF_8);
				streamWriter.write(entityWriter.toString());
				streamWriter.flush();
			}
		}
	}



	/**
	 * Instances of this class unmarshal JSON representations into generic objects. The resulting object will consist of
	 * primitive wrapper values (Boolean, Character, Long, Double only), strings, object arrays, and/or string-maps recursively
	 * containing the aforementioned types. Note that this class is designed to unmarshal raw JDBC table data from JSON, not
	 * JAX-B annotated entities (use MOXY instead).
	 */
	@Provider
	@Consumes(MediaType.APPLICATION_JSON)
	static private class GenericJsonEntityReader implements MessageBodyReader<Object> {

		/**
		 * {@inheritDoc}
		 */
		@Override
		public boolean isReadable (final Class<?> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType) {
			return MediaType.APPLICATION_JSON_TYPE.isCompatible(mediaType);
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public Object readFrom (final Class<Object> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType, final MultivaluedMap<String,String> httpHeaders, final InputStream entityStream) throws IllegalArgumentException, IOException, WebApplicationException {
			assert this.isReadable(type, genericType, annotations, mediaType);
			final BufferedReader reader = new BufferedReader(new InputStreamReader(entityStream, StandardCharsets.UTF_8), 1);
			return RestJsonMarshaller.unmarshal(reader);
		}
	}



	/**
	 * Instances of this class unmarshal JSON representations into generic lists. The resulting list will consist of primitive
	 * wrapper values (Boolean, Character, Long, Double only), strings, object arrays, and/or string-maps recursively containing
	 * the aforementioned types. Note that this class is designed to unmarshal raw JDBC table data from JSON, not JAX-B
	 * annotated entities (use MOXY instead).
	 */
	@Provider
	@Consumes(MediaType.APPLICATION_JSON)
	static private class GenericJsonListReader implements MessageBodyReader<List<Object>> {

		/**
		 * {@inheritDoc}
		 */
		@Override
		public boolean isReadable (final Class<?> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType) {
			return MediaType.APPLICATION_JSON_TYPE.isCompatible(mediaType);
		}


		/**
		 * {@inheritDoc}
		 */
		@Override
		public List<Object> readFrom (final Class<List<Object>> type, final Type genericType, final Annotation[] annotations, final MediaType mediaType, final MultivaluedMap<String,String> httpHeaders, final InputStream entityStream) throws IllegalArgumentException, IOException, WebApplicationException {
			assert this.isReadable(type, genericType, annotations, mediaType);
			final BufferedReader reader = new BufferedReader(new InputStreamReader(entityStream, StandardCharsets.UTF_8), 1);
			return RestJsonMarshaller.unmarshalCollection(reader);
		}
	}


	/**
	 * Marshals the given entity into a JSON representation and writes the result into the given writer.
	 * @param entity the entity
	 * @param entityWriter the entity writer
	 * @throws NullPointerException if the given entity writer is {@code null}
	 * @throws IOException if there is an I/O related problem
	 * @see Objects#toString(Object)
	 */
	@SuppressWarnings("unchecked")
	static private void marshal (final Object entity, final Writer entityWriter) throws NullPointerException, IOException {
		if (entity == null | entity instanceof Number | entity instanceof Boolean) {
			entityWriter.write(Objects.toString(entity));
		} else if (entity.getClass().isArray()) {
			marshalArray(entity, entityWriter);
		} else if (entity instanceof Collection) {
			marshalCollection((Collection<Object>) entity, entityWriter);
		} else if (entity instanceof Map) {
			marshalMap((Map<Object,Object>) entity, entityWriter);
		} else if (entity instanceof Character) {
			writeQuotedString(entity, '\'', entityWriter);
		} else {
			writeQuotedString(entity, '"', entityWriter);
		}
	}


	/**
	 * Marshals the given array into a JSON representation and writes the result into the given writer.
	 * @param entity the entity
	 * @param entityWriter the entity writer
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given entity is not an array
	 * @throws IOException if there is an I/O related problem
	 * @see Objects#toString(Object)
	 */
	static private void marshalArray (final Object entity, final Writer entityWriter) throws NullPointerException, IllegalArgumentException, IOException {
		if (!entity.getClass().isArray()) throw new IllegalArgumentException();

		entityWriter.write('[');
		for (int index = 0, stop = Array.getLength(entity); index < stop; ++index) {
			marshal(Array.get(entity, index), entityWriter);
			if (index != stop - 1) entityWriter.write(",");
		}
		entityWriter.write(']');
	}


	/**
	 * Marshals the given collection into a JSON representation and writes the result into the given writer.
	 * @param entity the entity
	 * @param entityWriter the entity writer
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IOException if there is an I/O related problem
	 * @see Objects#toString(Object)
	 */
	static private void marshalCollection (final Collection<?> entity, final Writer entityWriter) throws NullPointerException, IOException {
		entityWriter.write('[');
		for (final Iterator<?> iterator = entity.iterator(); iterator.hasNext();) {
			marshal(iterator.next(), entityWriter);
			if (iterator.hasNext()) entityWriter.write(",");
		}
		entityWriter.write(']');
	}


	/**
	 * Marshals the given string map into a JSON representation and writes the result into the given writer.
	 * @param entity the entity
	 * @param entityWriter the entity writer
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IOException if there is an I/O related problem
	 * @see Objects#toString(Object)
	 */
	static private void marshalMap (final Map<Object,Object> map, final Writer entityWriter) throws NullPointerException, IllegalArgumentException, IOException {
		entityWriter.write('{');

		for (final Iterator<Map.Entry<Object,Object>> iterator = map.entrySet().iterator(); iterator.hasNext();) {
			final Map.Entry<Object,Object> entry = iterator.next();
			writeQuotedString(entry.getKey(), '"', entityWriter);
			entityWriter.write(':');
			marshal(entry.getValue(), entityWriter);
			if (iterator.hasNext()) entityWriter.write(",");
		}

		entityWriter.write('}');
	}


	/**
	 * Writes the given entity's quoted string representation into the given writer.
	 * @param entity the entity
	 * @param quote the quote
	 * @param entityWriter the entity writer
	 * @throws NullPointerException if the given entity writer is {@code null}
	 * @throws IOException if there is an I/O related problem
	 * @see Objects#toString(Object)
	 */
	static private void writeQuotedString (final Object entity, final char quote, final Writer entityWriter) throws NullPointerException, IOException {
		entityWriter.write(quote);
		entityWriter.write(Objects.toString(entity));
		entityWriter.write(quote);
	}


	/**
	 * Unmarshals an object from the given reader, and returns it.
	 * @param entityReader the entity reader
	 * @return the unmarshaled object
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the reader's next characters do not represent a valid object
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws IOException if there is an I/O related problem, or if the given reader doesn't support mark/reset
	 */
	static private Object unmarshal (final Reader entityReader) throws NullPointerException, IllegalArgumentException, IOException {
		final int character = readNonWhitespaceCharacter(entityReader, true);
		switch (character) {
			case -1:
				throw new EOFException();
			case 'n':
				return unmarshalNull(entityReader);
			case 't':
			case 'f':
				return unmarshalBoolean(entityReader);
			case '\'':
				return unmarshalCharacter(entityReader);
			case '"':
				return unmarshalString(entityReader);
			case '[':
				return unmarshalCollection(entityReader);
			case '{':
				return unmarshalMap(entityReader);
			default:
				return unmarshalNumber(entityReader);
		}
	}


	/**
	 * Unmarshals a {@code null} value from the given reader, and returns it.
	 * @param entityReader the entity reader
	 * @return the unmarshaled null value
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws ProtocolException if the reader's next characters do not represent a valid null value
	 * @throws IOException if there is an I/O related problem
	 */
	static private Object unmarshalNull (final Reader entityReader) throws NullPointerException, IOException {
		final char[] characters = new char[4];
		for (int index = 0; index < characters.length; ++index) {
			final int character = entityReader.read();
			if (character == -1) throw new EOFException();
			characters[index] = (char) character;
		}

		if (!"null".equals(String.copyValueOf(characters))) throw new ProtocolException();
		return null;
	}


	/**
	 * Unmarshals a double-quoted string value from the given reader, and returns it (without quotes).
	 * @param entityReader the entity reader
	 * @return the unmarshaled String value
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the reader's next character is not a double-quote
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws IOException if there is an I/O related problem
	 */
	static private String unmarshalString (final Reader entityReader) throws NullPointerException, IllegalArgumentException, IOException {
		if (entityReader.read() != '"') throw new IllegalArgumentException();
		return readDelimitedString(entityReader, '"');
	}


	/**
	 * Unmarshals an unquoted boolean value from the given reader, and returns it.
	 * @param entityReader the entity reader
	 * @return the unmarshaled boolean value
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the reader's next characters do not represent a valid boolean value
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws IOException if there is an I/O related problem
	 */
	static private Boolean unmarshalBoolean (final Reader entityReader) throws NullPointerException, IllegalArgumentException, IOException {
		final String text = readDelimitedString(entityReader, 'e');
		switch (text) {
			case "true":
				return Boolean.TRUE;
			case "false":
				return Boolean.FALSE;
			default:
				throw new ProtocolException();
		}
	}


	/**
	 * Unmarshals a single-quoted character value from the given reader, and returns it (without quotes).
	 * @param entityReader the entity reader
	 * @return the unmarshaled character value
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws ProtocolException if the reader's next characters do not represent a valid character value
	 * @throws IOException if there is an I/O related problem
	 */
	static private Character unmarshalCharacter (final Reader entityReader) throws NullPointerException, IOException {
		if (entityReader.read() != '\'') throw new ProtocolException();
		final int character = entityReader.read();
		if (character == -1) throw new EOFException();
		if (entityReader.read() != '\'') throw new ProtocolException();

		return Character.valueOf((char) character);
	}


	/**
	 * Unmarshals an unquoted numeric value from the given reader, and returns it (without quotes).
	 * @param entityReader the entity reader
	 * @return the unmarshaled double or long value, depending on the presence of a decimal point
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the reader's next characters do not represent a valid double or long value
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws IOException if there is an I/O related problem, or if the given reader doesn't support mark/reset
	 */
	static private Number unmarshalNumber (final Reader entityReader) throws NullPointerException, IOException {
		try (StringWriter writer = new StringWriter()) {
			while (true) {
				entityReader.mark(1);
				final int character = entityReader.read();

				if ((character < '0' | character > '9') & character != '+' & character != '-' & character != '.' & character != 'e' & character != 'E') {
					entityReader.reset();
					final String text = writer.toString();
					return text.contains(".") ? (Number) Double.valueOf(text) : (Number) Long.valueOf(text);
				}

				writer.write(character);
			}
		}
	}


	/**
	 * Unmarshals a <tt>[]</tt>-bracketed list value from the given reader, and returns it (without brackets).
	 * @param entityReader the entity reader
	 * @return the unmarshaled Object list
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws ProtocolException if the reader's next characters do not represent a valid list
	 * @throws IOException if there is an I/O related problem, or if the given reader doesn't support mark/reset
	 */
	static private List<Object> unmarshalCollection (final Reader entityReader) throws NullPointerException, IOException {
		final List<Object> result = new ArrayList<>();
		while (true) {
			final int character = readNonWhitespaceCharacter(entityReader, false);
			switch (character) {
				case -1:
					throw new EOFException();
				case '[':
				case ',':
					result.add(unmarshal(entityReader));
					break;
				case ']':
					return result;
				default:
					throw new ProtocolException();
			}
		}
	}


	/**
	 * Unmarshals a <tt>{}</tt>-bracketed map value from the given reader, and returns it (without brackets).
	 * @param entityReader the entity reader
	 * @return the unmarshaled String-Object map
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws ProtocolException if the reader's next characters do not represent a valid map
	 * @throws IOException if there is an I/O related problem, or if the given reader doesn't support mark/reset
	 */
	static private Map<String,Object> unmarshalMap (final Reader entityReader) throws NullPointerException, IOException {
		final Map<String,Object> result = new HashMap<>();
		while (true) {
			final int character = readNonWhitespaceCharacter(entityReader, false);
			switch (character) {
				case -1:
					throw new EOFException();
				case '{':
				case ',':
					final String key;
					if (readNonWhitespaceCharacter(entityReader, true) == '"') {
						key = unmarshalString(entityReader);
						if (entityReader.read() != ':') throw new ProtocolException();
					} else {
						key = readDelimitedString(entityReader, ':').trim();
					}

					result.put(key, unmarshal(entityReader));
					break;
				case '}':
					return result;
				default:
					throw new ProtocolException();
			}
		}
	}


	/**
	 * Returns the given reader's next characters up until and excluding the given delimiter as a string, including whitespace.
	 * Note that this operation treats backslash-escaped characters as non-delimiters.
	 * @param entityReader the entity reader
	 * @param delimiter the text delimiter
	 * @return the string
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws EOFException if end-of-stream is reached prematurely
	 * @throws IOException if there is an I/O related problem
	 */
	static private String readDelimitedString (final Reader entityReader, final char delimiter) throws NullPointerException, IOException {
		try (StringWriter writer = new StringWriter()) {
			for (int predecessor = -1, successor = entityReader.read(); true; predecessor = successor, successor = entityReader.read()) {
				if (successor == -1) throw new EOFException();
				if (predecessor != '\\' & successor == delimiter) return writer.toString();
				writer.write(successor);
			}
		}
	}


	/**
	 * Skips any whitespace characters in the given reader, and returns it's next non-whitespace character, or -1 if the end of
	 * the reader has been reached.
	 * @param entityReader the entity reader
	 * @param peek whether or not to reset the reader after reading the next non-whitespace character
	 * @return the next non-whitespace character
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IOException if there is an I/O related problem, or if the given reader doesn't support mark/reset
	 */
	static private int readNonWhitespaceCharacter (final Reader entityReader, final boolean peek) throws NullPointerException, IOException {
		while (true) {
			if (peek) entityReader.mark(1);
			final int character = entityReader.read();
			if (character == -1 || !Character.isWhitespace((char) character)) {
				if (peek) entityReader.reset();
				return character;
			}
		}
	}
}