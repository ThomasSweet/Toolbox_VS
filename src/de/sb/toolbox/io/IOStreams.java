package de.sb.toolbox.io;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.io.StringWriter;
import java.io.Writer;
import java.net.SocketException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Stream;
import de.sb.toolbox.Copyright;


/**
 * This facade offers stream related convenience methods.Note that this class is declared final because it is a facade, and
 * therefore not supposed to be extended.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class IOStreams {

	/**
	 * Prevents instantiation.
	 */
	private IOStreams () {}


	/**
	 * Reads all remaining bytes from the given byte source, and writes them to the given byte sink. Returns the number of bytes
	 * copied, and closes neither source nor sink. Note that large copy buffers speed up processing, but consume more memory.
	 * Also note that {@link SocketException} is treated as a kind of EOF due to to other side terminating the stream.
	 * @param byteSource the byte source
	 * @param byteSink the byte sink
	 * @param bufferSize the buffer size, in number of bytes
	 * @return the number of bytes copied
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given buffer size is negative
	 * @throws IOException if there is an I/O related problem
	 */
	static public long copy (final InputStream byteSource, final OutputStream byteSink, final int bufferSize) throws IOException {
		if (bufferSize <= 0) throw new IllegalArgumentException();
		final byte[] buffer = new byte[bufferSize];

		long bytesCopied = 0;
		try {
			for (int bytesRead = byteSource.read(buffer); bytesRead != -1; bytesRead = byteSource.read(buffer)) {
				byteSink.write(buffer, 0, bytesRead);
				bytesCopied += bytesRead;
			}
		} catch (final SocketException exception) {
			// treat as EOF because a TCP stream has been closed by the other side
		}
		return bytesCopied;
	}


	/**
	 * Reads all remaining characters from the given char source, and writes them to the given char sink. Returns the number of
	 * characters copied, and closes neither source nor sink. Note that large copy buffers speed up processing, but consume more
	 * memory. Also note that {@link SocketException} is treated as a kind of EOF due to to other side terminating the stream.
	 * @param charSource the byte source
	 * @param charSink the byte sink
	 * @param bufferSize the buffer size, in number of characters
	 * @return the number of characters copied
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given buffer size is negative
	 * @throws IOException if there is an I/O related problem
	 */
	static public long copy (final Reader charSource, final Writer charSink, final int bufferSize) throws IOException {
		if (bufferSize <= 0) throw new IllegalArgumentException();
		final char[] buffer = new char[bufferSize];

		long charsCopied = 0;
		try {
			for (int charsRead = charSource.read(buffer); charsRead != -1; charsRead = charSource.read(buffer)) {
				charSink.write(buffer, 0, charsRead);
				charsCopied += charsRead;
			}
		} catch (final SocketException exception) {
			// treat as EOF because a TCP stream has been closed by the other side
		}
		return charsCopied;
	}


	/**
	 * Reads all remaining bytes from the given byte source, and returns them as a byte array.
	 * @param byteSource the byte source
	 * @return the bytes
	 * @throws NullPointerException if the given byte source is {@code null}
	 * @throws IOException if there is an I/O related problem
	 */
	static public byte[] read (final InputStream byteSource) throws IOException {
		try (ByteArrayOutputStream byteSink = new ByteArrayOutputStream()) {
			copy(byteSource, byteSink, 0x10000);
			return byteSink.toByteArray();
		}
	}


	/**
	 * Reads all remaining characters from the given char source, and returns them as a char array. If requested, the char
	 * source will be closed upon completion.
	 * @param charSource the char source
	 * @param closeUponCompletion whether or not to close the byte source upon completion
	 * @return the characters
	 * @throws NullPointerException if the given char source is {@code null}
	 * @throws IOException if there is an I/O related problem
	 */
	static public String read (final Reader charSource) throws IOException {
		try (StringWriter charSink = new StringWriter()) {
			copy(charSource, charSink, 0x10000);
			return charSink.toString();
		}
	}


	/**
	 * Returns the path and binary content of all the files within the given file system. Note that this operation is designed
	 * to work with virtual file systems.
	 * @param fileSystem the (virtual) file system
	 * @return the file names and their respective binary content as a map
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IOException if there is an I/O related problem
	 */
	static public SortedMap<String,byte[]> read (final FileSystem fileSystem) throws IOException {
		final SortedMap<String,byte[]> result = new TreeMap<>();

		for (final Path directory : fileSystem.getRootDirectories()) {
			try (Stream<Path> stream = Files.walk(directory)) {
				for (final Iterator<Path> iterator = stream.iterator(); iterator.hasNext();) {
					final Path path = iterator.next();
					if (Files.isRegularFile(path) & Files.isReadable(path)) result.put(path.toString(), Files.readAllBytes(path));
				}
			}
		}

		return result;
	}
}