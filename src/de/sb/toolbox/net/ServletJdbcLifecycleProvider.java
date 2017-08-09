package de.sb.toolbox.net;

import java.io.PrintWriter;
import java.sql.Connection;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.naming.InitialContext;
import javax.servlet.ServletContext;
import javax.servlet.ServletContextEvent;
import javax.servlet.ServletContextListener;
import javax.servlet.ServletRequestEvent;
import javax.servlet.ServletRequestListener;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;
import javax.sql.DataSource;
import de.sb.toolbox.Copyright;


/**
 * HTTP lifecycle provider managing request demarcated "physical" JDBC transactions. The class is designed to be configured as a
 * listener within a webapp deployment descriptor. Upon servlet context events, it creates/closes a JDBC data source. Upon
 * servlet request events, it begins/rolls back a JDBC transaction. The following context initialization parameters can be
 * specified within the deployment descriptor:
 * <ul>
 * <li>{@code de.sb.toolboxee.servlet.JDBC_DATA_SOURCE}: The data source JNDI name to be used for connection creation. This
 * parameter must be configured, there is no default.</li>
 * <li>{@code de.sb.toolboxee.servlet.JDBC_ISOLATION_LEVEL}: The transactional isolation level to be used for the connection.
 * Can be one of NONE, READ_UNCOMMITTED, READ_COMMITTED, REPEATABLE_READ, or SERIALIZABLE. Default is SERIALIZABLE.</li>
 * </ul>
 * <p>
 * The following request scoped attributes become available due to this listener:
 * </p>
 * <ul>
 * <li>{@code de.sb.toolboxee.servlet.JDBC_CONNECTION}: A JDBC connection that is initialized for transactions with isolation
 * level REPEATABLE_READ upon request begin. It is rolled back, closed, and removed upon request end.</li>
 * </ul>
 * Note that this implementation is based on servlet specification v2.5, SRV 10.3.3, stating that a container will only
 * instantiate a single instance of each configured listener per webapp context. Therefore, listener instance variables are used
 * to store configuration information.
 */
@Copyright(year = 2005, holders = "Sascha Baumeister")
public final class ServletJdbcLifecycleProvider implements ServletContextListener, ServletRequestListener {
	static private final String ATTRIBUTE_KEY_PREFIX = ServletJdbcLifecycleProvider.class.getPackage().getName() + ".";
	static private final String JDBC_DATA_SOURCE_KEY = ATTRIBUTE_KEY_PREFIX + "JDBC_DATA_SOURCE";
	static private final String JDBC_ISOLATION_LEVEL_KEY = ATTRIBUTE_KEY_PREFIX + "JDBC_ISOLATION_LEVEL";
	static public final String JDBC_CONNECTION_KEY = ATTRIBUTE_KEY_PREFIX + "JDBC_CONNECTION";



	static private enum IsolationLevel {
		NONE (Connection.TRANSACTION_NONE), READ_UNCOMMITTED (Connection.TRANSACTION_READ_UNCOMMITTED), READ_COMMITTED (Connection.TRANSACTION_READ_COMMITTED), REPEATABLE_READ (Connection.TRANSACTION_REPEATABLE_READ), SERIALIZABLE (Connection.TRANSACTION_SERIALIZABLE);

		private final int level;


		private IsolationLevel (final int level) {
			this.level = level;
		}
	}


	private DataSource dataSource = null;
	private IsolationLevel isolationLevel = null;


	/**
	 * Creates the receiver's JDBC data source and transactional isolation level from configuration data.
	 * @param event the application event
	 */
	public void contextInitialized (final ServletContextEvent event) {
		final ServletContext application = event.getServletContext();
		final String jdbcDataSourceName = application.getInitParameter(JDBC_DATA_SOURCE_KEY);
		try {
			this.dataSource = (DataSource) new InitialContext().lookup(jdbcDataSourceName);
			this.dataSource.setLogWriter(new PrintWriter(System.out));
			Logger.getGlobal().log(Level.CONFIG, "JDBC data source created: application={0}, jndiName={1}.", new Object[] { application.getContextPath(), jdbcDataSourceName });
		} catch (final Exception exception) {
			Logger.getGlobal().log(Level.SEVERE, "Unknown JDBC data source: application={0}, jndiName={1}.", new Object[] { application.getContextPath(), jdbcDataSourceName });
			throw new IllegalStateException(exception);
		}

		try {
			final String level = application.getInitParameter(JDBC_ISOLATION_LEVEL_KEY);
			this.isolationLevel = level == null ? null : IsolationLevel.valueOf(level);
			Logger.getGlobal().log(Level.CONFIG, "JDBC transactional isolation level set as configured: application={0}, isolationLevel={1}.", new Object[] { application.getContextPath(), this.isolationLevel });
		} catch (final NullPointerException exception) {
			this.isolationLevel = IsolationLevel.SERIALIZABLE;
			Logger.getGlobal().log(Level.CONFIG, "JDBC transactional isolation level set from default: application={0}, isolationLevel={1}.", new Object[] { application.getContextPath(), this.isolationLevel });
		} catch (final IllegalArgumentException exception) {
			this.isolationLevel = IsolationLevel.SERIALIZABLE;
			Logger.getGlobal().log(Level.WARNING, "Illegal JDBC transactional isolation level configured, set to default: application={0}, isolationLevel={1}.", new Object[] { application.getContextPath(), this.isolationLevel });
		}
	}


	/**
	 * Releases the receiver's data source.
	 * @param event the application event
	 */
	public void contextDestroyed (final ServletContextEvent event) {
		this.dataSource = null;
	}


	/**
	 * Starts a new JDBC request demarcated transaction.
	 * @param event the request event
	 */
	public void requestInitialized (final ServletRequestEvent event) {
		if (event.getServletRequest() instanceof HttpServletRequest) {
			final HttpServletRequest request = (HttpServletRequest) event.getServletRequest();
			final HttpSession session = request.getSession();
			try {
				final Connection connection = this.dataSource.getConnection();
				connection.setAutoCommit(false);
				if (this.isolationLevel != null) connection.setTransactionIsolation(this.isolationLevel.level);
				request.setAttribute(JDBC_CONNECTION_KEY, connection);
			} catch (final Exception exception) {
				Logger.getGlobal().log(Level.SEVERE, "Cannot create JDBC connection: application={0}, session={1}.", new Object[] { session.getServletContext().getContextPath(), session.getId() });
				throw new IllegalStateException(exception);
			}
		}
	}


	/**
	 * Rolls back the JDBC request demarcated transaction.
	 * @param event the request event
	 */
	public void requestDestroyed (final ServletRequestEvent event) {
		if (event.getServletRequest() instanceof HttpServletRequest) {
			final HttpServletRequest request = (HttpServletRequest) event.getServletRequest();
			final HttpSession session = request.getSession();
			try {
				final Connection connection = (Connection) request.getAttribute(JDBC_CONNECTION_KEY);
				connection.rollback();
				connection.close();
			} catch (final Exception exception) {
				Logger.getGlobal().log(Level.WARNING, "Cannot close JDBC connection: application={0}, session={1}.", new Object[] { session.getServletContext().getContextPath(), session.getId() });
			}
		}
	}
}