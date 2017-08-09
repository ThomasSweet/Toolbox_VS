package de.sb.toolbox.io;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.channels.AsynchronousCloseException;
import java.nio.channels.ClosedChannelException;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.FileLockInterruptionException;
import java.nio.channels.InterruptibleChannel;
import java.nio.channels.NonReadableChannelException;
import java.nio.channels.NonWritableChannelException;
import java.nio.channels.OverlappingFileLockException;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.FileAttribute;
import java.util.Set;
import de.sb.toolbox.Copyright;


/**
 * Instances of this class provide combined read&write access to an underlying seekable resource, usually a file. They provide
 * two stream views for read and write access, which both operate on the same position pointer. Therefore, reading a couple of
 * bytes has the effect the the current write pointer also moves forward, and vice versa. This implies that the input stream
 * view is used for both reading data, and repositioning the underlying channel, while the output stream view is used solely for
 * writing.<br />
 * Note the main advantage of this class over directly accessing a seekable byte channel is the avoidance of having to deal with
 * the {@link ByteBuffer} API. Being forced to deal with the latter in every read and write operation is probably the main
 * reason why the whole channel API is less widely known and used than it deserves, likely because the ByteBuffer-API is more
 * suited for data conversion, assemblage, and buffering than streaming I/O.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class StreamChannel implements AutoCloseable, InterruptibleChannel {

	private final FileChannel fileChannel;
	private final FileChannelInputStream fileChannelSource;
	private final FileChannelOutputStream fileChannelSink;


	/**
	 * Creates a new stream channel for the given path and open options. Note that {@link StandardOpenOption#APPEND APPEND} is
	 * not allowed because is seriously messes with positioning in conjunction with write methods.
	 * @param path the path
	 * @param options the open options
	 * @return the new stream channel
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if one of the given options is {@link StandardOpenOption#APPEND APPEND}
	 * @throws IOException if there is an I/O related problem
	 */
	static public StreamChannel newStreamChannel (final Path path, final OpenOption... options) throws IOException {
		for (final OpenOption option : options) {
			if (option.equals(StandardOpenOption.APPEND)) throw new IllegalArgumentException();
		}

		final FileChannel fileChannel = FileChannel.open(path, options);
		return new StreamChannel(fileChannel);
	}


	/**
	 * Creates a new stream channel for the given path, open options and file attributes. Note that
	 * {@link StandardOpenOption#APPEND APPEND} is not allowed because is seriously messes with positioning in conjunction with
	 * write methods.
	 * @param path the path
	 * @param options the open options
	 * @param attributes the desired file attributes
	 * @return the new stream channel
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if one of the given options is {@link StandardOpenOption#APPEND APPEND}
	 * @throws IOException if there is an I/O related problem
	 */
	static public StreamChannel newStreamChannel (final Path path, final Set<? extends OpenOption> options, final FileAttribute<?>... attributes) throws IOException {
		for (final OpenOption option : options) {
			if (option.equals(StandardOpenOption.APPEND)) throw new IllegalArgumentException();
		}

		final FileChannel fileChannel = FileChannel.open(path, options, attributes);
		return new StreamChannel(fileChannel);
	}


	/**
	 * Creates a new instance based on the given file channel
	 * @param fileChannel the file channel
	 * @throws IOException if there is an I/O problem forcing updates to
	 * @throws NullPointerException if the given file channel is {@code null}
	 */
	protected StreamChannel (final FileChannel fileChannel) throws IOException {
		if (fileChannel == null) throw new NullPointerException();

		this.fileChannel = fileChannel;
		this.fileChannelSource = new FileChannelInputStream();
		this.fileChannelSink = new FileChannelOutputStream();
	}


	/**
	 * Forces any updates of the underlying channel's file to be written to the storage device that contains it. Then closes the
	 * underlying channel, and therefore also both stream views. Does nothing if this channel is already closed, or in the
	 * process of closing.
	 * @throws IOException if there is an I/O related problem
	 */
	public void close () throws IOException {
		try {
			StreamChannel.this.fileChannel.force(false);
		} catch (final Exception exception) {}
		this.fileChannel.close();
	}


	/**
	 * Returns {@code true} if the underlying channel is open, {@code false} if closed.
	 * @return whether or not the underlying channel is open
	 */
	public boolean isOpen () {
		return this.fileChannel.isOpen();
	}


	/**
	 * Returns the size of the underlying channel.
	 * @return the size
	 * @throws ClosedChannelException If this channel is closed
	 * @throws IOException if an I/O related problem occurs
	 */
	public long getSize () throws IOException {
		return this.fileChannel.size();
	}


	/**
	 * Sets the size of the underlying channel without affecting the underlying channel's position. If the channel size is
	 * shrunken, this position may become beyond EOF.
	 * @param size the size
	 * @throws IllegalArgumentException if the given size is strictly negative
	 * @throws NonWritableChannelException if the underlying channel was not opened for writing
	 * @throws ClosedChannelException If this channel is closed
	 * @throws IOException if an I/O related problem occurs
	 */
	public void setSize (final long size) throws IOException {
		synchronized (StreamChannel.this.fileChannel) {
			final long position = this.fileChannel.position();
			try {
				if (this.fileChannel.size() >= size) {
					this.fileChannel.truncate(size);
				} else {
					this.fileChannel.write(ByteBuffer.wrap(new byte[1]), size - 1);
				}
			} finally {
				this.fileChannel.position(position);
			}
		}
	}


	/**
	 * Acquires an exclusive lock on the underlying channel's file. An invocation of this method will block until the whole file
	 * can be locked, this channel is closed, or the invoking thread is interrupted, whichever comes first. File locks are held
	 * on behalf of the entire Java virtual machine. They are not suitable for controlling access to a file by multiple threads
	 * within the same virtual machine.
	 * @return a lock object representing the newly-acquired lock
	 * @throws NonWritableChannelException if the underlying channel was not opened for writing
	 * @throws FileLockInterruptionException if the invoking thread is interrupted while waiting to acquire the lock, or if the
	 *         current thread's interrupt status is already set when this method is invoked
	 * @throws ClosedChannelException if this channel is closed
	 * @throws AsynchronousCloseException if another thread closes this channel while the invoking thread is blocked in this
	 *         method
	 * @throws OverlappingFileLockException if a lock on any file region is already held by this Java virtual machine, or if
	 *         another thread is already blocked in this method and is attempting to lock any region of the same file
	 * @throws IOException if some other I/O error occurs
	 * @see FileChannel#lock()
	 * @see FileLock#release()
	 */
	public FileLock lock () throws IOException {
		return this.fileChannel.lock();
	}


	/**
	 * Acquires a lock on the given region of the underlying channel's file. An invocation of this method will block until the
	 * region can be locked, this channel is closed, or the invoking thread is interrupted, whichever comes first. File locks
	 * are held on behalf of the entire Java virtual machine. They are not suitable for controlling access to a file by multiple
	 * threads within the same virtual machine.
	 * @param position The position at which the locked region is to start; must be non-negative
	 * @param size the size of the locked region; must be non-negative, and the sum position + size must be non-negative
	 * @param shared {@code true} to request a shared lock, {@code false} to request an exclusive lock
	 * @return a lock object representing the newly-acquired lock
	 * @throws IllegalArgumentException if the preconditions on the parameters do not hold
	 * @throws NonReadableChannelException if shared is {@code true} but the underlying channel was not opened for reading
	 * @throws NonWritableChannelException if shared is {@code false} but the underlying channel was not opened for writing
	 * @throws FileLockInterruptionException if the invoking thread is interrupted while waiting to acquire the lock, or if the
	 *         current thread's interrupt status is already set when this method is invoked
	 * @throws ClosedChannelException if this channel is closed
	 * @throws AsynchronousCloseException if another thread closes this channel while the invoking thread is blocked in this
	 *         method
	 * @throws OverlappingFileLockException if a lock that overlaps the requested region is already held by this Java virtual
	 *         machine, or if another thread is already blocked in this method and is attempting to lock an overlapping region
	 *         of the same file
	 * @throws IOException if some other I/O error occurs
	 * @see FileChannel#lock(long, long, boolean)
	 * @see FileLock#release()
	 */
	public FileLock lock (final long position, final long size, final boolean shared) throws IOException {
		return this.fileChannel.lock(position, size, shared);
	}


	/**
	 * Attempts to acquire an exclusive lock on the underlying channel's file. This method does not block. An invocation always
	 * returns immediately, either having acquired the exclusive lock on the file or having failed to do so. If it fails to
	 * acquire a lock because a lock is already held by another program then it returns {@code null}. If it fails to acquire a
	 * lock for any other reason then an appropriate exception is thrown. File locks are held on behalf of the entire Java
	 * virtual machine. They are not suitable for controlling access to a file by multiple threads within the same virtual
	 * machine.
	 * @return a lock object representing the newly-acquired lock, or {@code null} if the lock could not be acquired because
	 *         another program holds an overlapping lock
	 * @throws NonWritableChannelException if the underlying channel was not opened for writing
	 * @throws ClosedChannelException if this channel is closed
	 * @throws OverlappingFileLockException if a lock is already held by this Java virtual machine, or if another thread is
	 *         already blocked in this method and is attempting to lock some region of the same file
	 * @throws IOException if some other I/O error occurs
	 * @see FileChannel#tryLock()
	 * @see FileLock#release()
	 */
	public FileLock tryLock () throws IOException {
		return this.fileChannel.tryLock();
	}


	/**
	 * Attempts to acquire a lock on the given region of the underlying channel's file. This method does not block. An
	 * invocation always returns immediately, either having acquired a lock on the requested region or having failed to do so.
	 * If it fails to acquire a lock because an overlapping lock is held by another program then it returns null. If it fails to
	 * acquire a lock for any other reason then an appropriate exception is thrown. File locks are held on behalf of the entire
	 * Java virtual machine. They are not suitable for controlling access to a file by multiple threads within the same virtual
	 * machine.
	 * @param position The position at which the locked region is to start; must be non-negative
	 * @param size the size of the locked region; must be non-negative, and the sum position + size must be non-negative
	 * @param shared {@code true} to request a shared lock, {@code false} to request an exclusive lock
	 * @return a lock object representing the newly-acquired lock, or {@code null} if the lock could not be acquired because
	 *         another program holds an overlapping lock
	 * @throws IllegalArgumentException if the preconditions on the parameters do not hold
	 * @throws NonReadableChannelException if shared is {@code true} but the underlying channel was not opened for reading
	 * @throws NonWritableChannelException if shared is {@code false} but the underlying channel was not opened for writing
	 * @throws ClosedChannelException if this channel is closed
	 * @throws OverlappingFileLockException if a lock that overlaps the requested region is already held by this Java virtual
	 *         machine, or if another thread is already blocked in this method and is attempting to lock an overlapping region
	 *         of the same file
	 * @throws IOException if some other I/O error occurs
	 * @see FileChannel#tryLock(long, long, boolean)
	 * @see FileLock#release()
	 */
	public FileLock tryLock (long position, long size, boolean shared) throws IOException {
		return this.fileChannel.tryLock(position, size, shared);
	}


	/**
	 * Returns the main input stream view of this channel. Repeated calls for this view will always return the same object. The
	 * resulting view (as well as it's writing counterpart) relies on this channel's positional state; both views are therefore
	 * designed to work as a tandem, usually within a single thread. This view is responsible for reading data, and
	 * fast-forwarding/rewinding the underlying channel's position; it's counterpart is responsible for writing data, and
	 * enforcing the persistence of changes. It doesn't necessarily need not be closed; if any channel view is closed, all other
	 * views are closed along with the associated channel.<br />
	 * In addition to basic input stream behavior, instances of this view will throw a {@link NonReadableChannelException}
	 * whenever asked to perform read I/O while the underlying channel was opened without the {@link StandardOpenOption#READ
	 * READ} option.<br />
	 * Note that any stream I/O operations except {@link InputStream#markSupported() markSupported()},
	 * {@link InputStream#mark(int) mark(int)}, and {@link InputStream#close() close()} will fail with a
	 * {@link ClosedChannelException} once the underlying channel is closed, or is in the process of closing. Also note that the
	 * {@link InputStream#skip(long) skip(long)} operation will always skip the given amount of bytes (unless it is negative),
	 * and may skip beyond the underlying channel's EOF; this is helpful when extending a file by a subsequent write I/O
	 * operation.
	 * @return this channel's input stream view
	 * @see #newInputStream(long)
	 * @see #getOutputStream()
	 */
	public InputStream getInputStream () {
		return this.fileChannelSource;
	}


	/**
	 * Returns the main output stream view of this channel. Repeated calls for this view will always return the same object. The
	 * resulting view (as well as it's reading counterpart) relies on this channel's positional state; both views are therefore
	 * designed to work as a tandem, usually within a single thread. This view is responsible for writing data, and enforcing
	 * the persistence of changes; it's counterpart is responsible for reading data, and fast- forwarding/rewinding the
	 * underlying channel's position. It doesn't necessarily need not be closed; if any channel view is closed, all other views
	 * are closed along with the associated channel.<br />
	 * In addition to basic output stream behavior, instances of this view will throw a {@link NonWritableChannelException}
	 * whenever asked to perform write I/O while the underlying channel was opened without the {@link StandardOpenOption#WRITE
	 * WRITE} option.<br />
	 * Note that any stream I/O operations except {@link OutputStream#close() close()} will fail with a
	 * {@link ClosedChannelException} once the underlying channel is closed, or is in the process of closing. Also note that the
	 * {@link OutputStream#flush() flush()} operation will force the underlying channel's updates to be written into the
	 * associated data store.
	 * @return this channel's output stream view
	 * @see #newOutputStream(long)
	 * @see #getInputStream()
	 */
	public OutputStream getOutputStream () {
		return this.fileChannelSink;
	}


	/**
	 * Returns a new input stream that is a detached view of this channel, owning a separate position counter. These views are
	 * primarily intended for multi-threaded read operations, with one thread per view. They dont't necessarily need not be
	 * closed; if any channel view is closed, all other views are closed along with the associated channel.<br />
	 * In addition to basic input stream behavior, instances of this view will throw a {@link NonReadableChannelException}
	 * whenever asked to perform read I/O while the underlying channel was opened without the {@link StandardOpenOption#READ
	 * READ} option.<br />
	 * Note that any stream I/O operations except {@link InputStream#markSupported() markSupported()},
	 * {@link InputStream#mark(int) mark(int)}, and {@link InputStream#close() close()} will fail with a
	 * {@link ClosedChannelException} once the underlying channel is closed, or is in the process of closing. Also note that the
	 * {@link InputStream#skip(long) skip(long)} operation will always skip the given amount of bytes (unless it is negative),
	 * and may skip beyond the underlying channel's EOF.
	 * @param position the initial position of the new view
	 * @return the detached channel input stream view created
	 * @throws IllegalArgumentException if the given position is negative
	 * @see #getInputStream()
	 * @see #newOutputStream(long)
	 */
	public InputStream newInputStream (final long position) {
		return new DetachedFileChannelInputStream(position);
	}


	/**
	 * Returns a new output stream that is a detached view of this channel, owning a separate position counter. These views are
	 * primarily intended for multi-threaded write operations, with one thread per view. They dont't necessarily need not be
	 * closed; if any channel view is closed, all other views are closed along with the associated channel.<br />
	 * In addition to basic output stream behavior, instances of this view will throw a {@link NonWritableChannelException}
	 * whenever asked to perform write I/O while the underlying channel was opened without the {@link StandardOpenOption#WRITE
	 * WRITE} option.<br />
	 * Note that any stream I/O operations except {@link OutputStream#close() close()} will fail with a
	 * {@link ClosedChannelException} once the underlying channel is closed, or is in the process of closing. Also note that the
	 * {@link OutputStream#flush() flush()} operation will force the underlying channel's updates to be written into the
	 * associated data store.
	 * @param position the initial position of the new view
	 * @return the detached channel output stream view created
	 * @throws IllegalArgumentException if the given position is negative
	 * @see #getOutputStream()
	 * @see #newInputStream(long)
	 */
	public OutputStream newOutputStream (final long position) {
		return new DetachedFileChannelOutputStream(position);
	}



	/**
	 * A channel's main input stream view for reading and repositioning, intended for single-threaded operations on a channel,
	 * in concert with it's counterpart that is responsible for write operations and enforcing persistence. In addition to basic
	 * input stream behavior, instances of this view will throw a {@link NonReadableChannelException} whenever asked to perform
	 * read I/O while the underlying channel was opened without the necessary {@link StandardOpenOption#READ READ} option.<br />
	 * Note that any stream I/O operations except {@link #markSupported()}, {@link #mark(int)}, and {@link #close()} will fail
	 * with a {@link ClosedChannelException} once the underlying channel is closed, or is in the process of closing. Also note
	 * that the {@link #skip(long)} operation will always skip the given amount of bytes (unless it is negative), and may skip
	 * beyond the underlying channel's {@code EOF}; this is helpful when extending a file by a subsequent write I/O operation.
	 */
	private class FileChannelInputStream extends InputStream {
		private volatile long markPosition = 0;


		@Override
		public int read () throws IOException {
			final byte[] buffer = new byte[1];
			final int bytesRead = this.read(buffer, 0, buffer.length);
			return bytesRead < 0 ? bytesRead : buffer[0] & 0xFF;
		}


		@Override
		public int read (final byte[] buffer, final int offset, final int length) throws IOException {
			if (length == 0) return 0;

			final ByteBuffer byteBuffer = ByteBuffer.wrap(buffer, offset, length);
			return StreamChannel.this.fileChannel.read(byteBuffer);
		}


		@Override
		public long skip (long offset) throws IOException {
			if (offset <= 0) return 0;

			synchronized (StreamChannel.this.fileChannel) {
				final long position = StreamChannel.this.fileChannel.position();
				StreamChannel.this.fileChannel.position(position + offset);
			}
			return offset;
		}


		@Override
		public int available () throws IOException {
			final long size = StreamChannel.this.fileChannel.size();
			final long position = StreamChannel.this.fileChannel.position();
			return (int) Math.min(Integer.MAX_VALUE, Math.max(0, size - position));
		}


		@Override
		public boolean markSupported () {
			return true;
		}


		@Override
		public void mark (final int readLimit) {
			try {
				this.markPosition = StreamChannel.this.fileChannel.position();
			} catch (final ClosedChannelException exception) {
				// do nothing, as any further operations will fail
			} catch (final IOException exception) {
				// try to make sure any further operations will fail
				try {
					StreamChannel.this.close();
				} catch (final Exception nestedException) {}
			}
		}


		@Override
		public void reset () throws IOException {
			StreamChannel.this.fileChannel.position(this.markPosition);
		}


		@Override
		public void close () throws IOException {
			StreamChannel.this.close();
		}
	}



	/**
	 * A channel's main output stream view for writing and enforcing persistence, intended for single-threaded operations on a
	 * channel, in concert with it's counterpart that is responsible for read and repositioning operations. In addition to basic
	 * output stream behavior, instances of this view will throw a {@link NonWritableChannelException} whenever asked to perform
	 * write I/O while the underlying channel was opened without the {@link StandardOpenOption#WRITE WRITE} option.<br />
	 * Note that any stream I/O operations except {@link #close()} will fail with a {@link ClosedChannelException} once the
	 * underlying channel is closed, or is in the process of closing. Also note that the {@link #flush()} operation will force
	 * the underlying channel's updates to be written into the associated data store.
	 */
	private class FileChannelOutputStream extends OutputStream {

		@Override
		public void write (final int unsignedByte) throws IOException {
			final byte[] buffer = new byte[1];
			buffer[0] = (byte) unsignedByte;
			this.write(buffer, 0, buffer.length);
		}


		@Override
		public void write (final byte[] buffer, final int offset, final int length) throws IOException {
			final ByteBuffer byteBuffer = ByteBuffer.wrap(buffer, offset, length);
			StreamChannel.this.fileChannel.write(byteBuffer);
		}


		@Override
		public void flush () throws IOException {
			StreamChannel.this.fileChannel.force(false);
		}


		@Override
		public void close () throws IOException {
			StreamChannel.this.close();
		}
	}



	/**
	 * A channel's detached input stream view for read operations, based on it's own position counter, and intended for
	 * multi-threaded read operations on a channel. In addition to basic input stream behavior, instances of this view will
	 * throw a {@link NonReadableChannelException} whenever asked to perform read I/O while the underlying channel was opened
	 * without the necessary {@link StandardOpenOption#READ READ} option.<br />
	 * Note that any stream I/O operations except {@link #markSupported()}, {@link #mark(int)}, and {@link #close()} will fail
	 * with a {@link ClosedChannelException} once the underlying channel is closed, or is in the process of closing. Also note
	 * that the {@link #skip(long)} operation will always skip the given amount of bytes (unless it is negative), and may skip
	 * beyond the underlying channel's {@code EOF}.
	 */
	private class DetachedFileChannelInputStream extends InputStream {
		private final Object mutex;
		private volatile long position;
		private volatile long markPosition;


		public DetachedFileChannelInputStream (final long position) {

			if (this.position < 0) throw new IllegalArgumentException();

			this.mutex = new Object();
			this.position = position;
			this.markPosition = 0;
		}


		@Override
		public int read () throws IOException {
			final byte[] buffer = new byte[1];
			final int bytesRead = this.read(buffer, 0, buffer.length);
			return bytesRead < 0 ? bytesRead : buffer[0] & 0xFF;
		}


		@Override
		public int read (final byte[] buffer, final int offset, final int length) throws IOException {
			if (length == 0) return 0;

			final ByteBuffer byteBuffer = ByteBuffer.wrap(buffer, offset, length);
			synchronized (this.mutex) {
				final int bytesRead = StreamChannel.this.fileChannel.read(byteBuffer, this.position);
				if (bytesRead != -1) this.position += bytesRead;
				return bytesRead;
			}
		}


		@Override
		public long skip (long offset) {
			if (offset <= 0) return 0;
			synchronized (this.mutex) {
				this.position += offset;
			}
			return offset;
		}


		@Override
		public int available () throws IOException {
			final long size = StreamChannel.this.fileChannel.size();
			return (int) Math.min(Integer.MAX_VALUE, Math.max(0, size - this.position));
		}


		@Override
		public boolean markSupported () {
			return true;
		}


		@Override
		public void mark (final int readLimit) {
			this.markPosition = this.position;
		}


		@Override
		public void reset () {
			this.position = this.markPosition;
		}


		@Override
		public void close () throws IOException {
			StreamChannel.this.close();
		}
	}



	/**
	 * A channel's detached output stream view for writing and force operations, based on it's own position counter, and
	 * intended for multi-threaded write operations on a channel. In addition to basic output stream behavior, instances of this
	 * view will throw a {@link NonWritableChannelException} whenever asked to perform write I/O while the underlying channel
	 * was opened without the {@link StandardOpenOption#WRITE WRITE} option.<br />
	 * Note that any stream I/O operations except {@link #close()} will fail with a {@link ClosedChannelException} once the
	 * underlying channel is closed, or is in the process of closing. Also note that the {@link #flush()} operation will force
	 * the underlying channel's updates to be written into the associated data store.
	 */
	private class DetachedFileChannelOutputStream extends OutputStream {
		private final Object mutex;
		private volatile long position;


		public DetachedFileChannelOutputStream (final long position) {

			if (this.position < 0) throw new IllegalArgumentException();

			this.mutex = new Object();
			this.position = position;
		}


		@Override
		public void write (final int unsignedByte) throws IOException {
			final byte[] buffer = new byte[1];
			buffer[0] = (byte) unsignedByte;
			this.write(buffer, 0, buffer.length);
		}


		@Override
		public void write (final byte[] buffer, final int offset, final int length) throws IOException {
			final ByteBuffer byteBuffer = ByteBuffer.wrap(buffer, offset, length);
			synchronized (this.mutex) {
				StreamChannel.this.fileChannel.write(byteBuffer, this.position);
				this.position += length;
			}
		}


		@Override
		public void flush () throws IOException {
			StreamChannel.this.fileChannel.force(false);
		}


		@Override
		public void close () throws IOException {
			StreamChannel.this.close();
		}
	}
}