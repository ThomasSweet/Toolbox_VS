package de.sb.toolbox.io;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import de.sb.toolbox.Copyright;


/**
 * Simple pipe-based socket-like objects that offer TCP-like two-way intra-process communications without requiring any
 * resources from the TCP/IP stack.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class PipeSocket implements Closeable {

	private final PipedInputStream inputStream;
	private final PipedOutputStream outputStream;
	private volatile boolean closed;


	/**
	 * Returns a pair of interconnected pipe sockets for intraprocess communication.
	 * @return the pair of pipe sockets
	 * @throws IOException if there is an I/O related problem
	 */
	static public PipeSocket[] newConnection () throws IOException {
		final PipedInputStream leftByteSource = new PipedInputStream();
		final PipedInputStream rightByteSource = new PipedInputStream();
		final PipedOutputStream leftByteSink = new PipedOutputStream(rightByteSource);
		final PipedOutputStream rightByteSink = new PipedOutputStream(leftByteSource);

		return new PipeSocket[] { new PipeSocket(leftByteSource, leftByteSink), new PipeSocket(rightByteSource, rightByteSink)
		};
	}


	/**
	 * Creates a new instance from a pair of piped streams.
	 * @param inputStream the piped input stream
	 * @param outputStream the piped output stream
	 * @throws NullPointerException if any of the given streams is {@code null}
	 */
	protected PipeSocket (final PipedInputStream inputStream, final PipedOutputStream outputStream) {
		if (inputStream == null | outputStream == null) throw new NullPointerException();
		this.inputStream = inputStream;
		this.outputStream = outputStream;
		this.closed = false;
	}


	/**
	 * {@inheritDoc}
	 * @throws IOException {@inheritDoc}
	 */
	@Override
	public void close () throws IOException {
		if (!this.closed) {
			this.closed = true;

			try {
				this.inputStream.close();
			} finally {
				this.outputStream.close();
			}
		}
	}


	/**
	 * Return the input stream
	 * @return the input stream
	 */
	public InputStream getInputStream () {
		return this.inputStream;
	}


	/**
	 * Return the output stream
	 * @return the output stream
	 */
	public OutputStream getOutputStream () {
		return this.outputStream;
	}


	/**
	 * Returns {@code true} if this socket is closed, {@code false} otherwise.
	 * @return whether or not this socket is closed
	 */
	public boolean isClosed () {
		return this.closed;
	}
}