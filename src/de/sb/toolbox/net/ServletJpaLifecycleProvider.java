package de.sb.toolbox.net;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.persistence.EntityManager;
import javax.persistence.EntityManagerFactory;
import javax.persistence.Persistence;
import javax.servlet.ServletContext;
import javax.servlet.ServletContextEvent;
import javax.servlet.ServletContextListener;
import javax.servlet.ServletRequestEvent;
import javax.servlet.ServletRequestListener;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;
import javax.servlet.http.HttpSessionEvent;
import javax.servlet.http.HttpSessionListener;
import de.sb.toolbox.Copyright;


/**
 * HTTP lifecycle provider managing session demarcated JPA "logical" business transactions based on request demarcated
 * "physical" JDBC transactions. The class is designed to be configured as a listener within a webapp deployment descriptor.
 * Upon servlet context events, it creates/closes an entity manager factory. Upon HTTP session events, it creates/closes an
 * entity manager that is published within session scope. Upon servlet request events, it begins/rolls back a JPA transaction.
 * The following context initialization parameters have to be specified within the deployment descriptor:
 * <ul>
 * <li>{@code de.sb.toolboxee.servlet.JDBC_DATA_SOURCE}: The data source JNDI name to be used for connection creation. There is
 * no default.</li>
 * <li>{@code de.sb.toolboxee.servlet.JPA_PERSISTENCE_UNIT}: The persistence unit name to be used for entity managers. There is
 * no default.</li>
 * </ul>
 * <p>
 * The following session scoped attributes become available due to this listener:
 * </p>
 * <ul>
 * <li>{@code de.sb.toolboxee.servlet.JPA_ENTITY_MANAGER}: A JPA entity manager primed with the configured (non JTA) data
 * source.</li>
 * </ul>
 * Note that this implementation is based on servlet specification v2.5, SRV 10.3.3, stating that a container will only
 * instantiate a single instance of each configured listener per webapp context. Therefore, listener instance variables are used
 * to store configuration information.
 */
@Copyright(year = 2005, holders = "Sascha Baumeister")
public final class ServletJpaLifecycleProvider implements ServletContextListener, HttpSessionListener, ServletRequestListener {
	static private final String ATTRIBUTE_KEY_PREFIX = ServletJpaLifecycleProvider.class.getPackage().getName() + ".";
	static private final String JDBC_DATA_SOURCE_KEY = ATTRIBUTE_KEY_PREFIX + "JDBC_DATA_SOURCE";
	static private final String JPA_PERSISTENCE_UNIT_KEY = ATTRIBUTE_KEY_PREFIX + "JPA_PERSISTENCE_UNIT";
	static public final String JPA_ENTITY_MANAGER_KEY = ATTRIBUTE_KEY_PREFIX + "JPA_ENTITY_MANAGER";

	private EntityManagerFactory entityManagerFactory = null;


	/**
	 * Creates the receiver's JPA entity manager factory from configuration data.
	 * @param event the application event
	 */
	public void contextInitialized (final ServletContextEvent event) {
		final ServletContext application = event.getServletContext();
		final String jdbcDataSourceName = application.getInitParameter(JDBC_DATA_SOURCE_KEY);
		final String jpaPersistenceUnitName = application.getInitParameter(JPA_PERSISTENCE_UNIT_KEY);

		final Map<String,String> configOverrides = new HashMap<String,String>();
		configOverrides.put("javax.persistence.nonJtaDataSource", jdbcDataSourceName);
		configOverrides.put("hibernate.show_sql", Boolean.toString(Logger.getGlobal().isLoggable(Level.FINER)));
		configOverrides.put("hibernate.use_sql_comments", Boolean.toString(Logger.getGlobal().isLoggable(Level.FINEST)));
		try {
			this.entityManagerFactory = Persistence.createEntityManagerFactory(jpaPersistenceUnitName, configOverrides);
			Logger.getGlobal().log(Level.CONFIG, "JPA entity manager factory created: application={0}, dataSourceJndiName={1}.", new Object[] { application.getContextPath(), jdbcDataSourceName });
		} catch (final Exception exception) {
			Logger.getGlobal().log(Level.SEVERE, "Creation of entity manager factory failed: application={0}, dataSourceJndiName={1}, message={2}.", new Object[] { application.getContextPath(), jdbcDataSourceName, exception.getMessage() });
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * Closes the receiver's JPA entity manager factory.
	 * @param event the application event
	 */
	public void contextDestroyed (final ServletContextEvent event) {
		this.entityManagerFactory.close();
		this.entityManagerFactory = null;
	}


	/**
	 * Creates a new session scoped JPA entity manager from the receiver's JPA entity manager factory.
	 * @param event the session event
	 */
	public void sessionCreated (final HttpSessionEvent event) {
		final HttpSession session = event.getSession();
		try {
			session.setAttribute(JPA_ENTITY_MANAGER_KEY, this.entityManagerFactory.createEntityManager());
		} catch (final Exception exception) {
			Logger.getGlobal().log(Level.SEVERE, "Creation of entity manager failed: application={0}, session={1}, factory={2}.", new Object[] { session.getServletContext().getContextPath(), session.getId(), this.entityManagerFactory });
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * Rolls back any remaining transaction and closes the session scoped JPA entity manager.
	 * @param event the session event
	 */
	public void sessionDestroyed (final HttpSessionEvent event) {
		final HttpSession session = event.getSession();
		final EntityManager entityManager = (EntityManager) session.getAttribute(JPA_ENTITY_MANAGER_KEY);
		try {
			if (entityManager.getTransaction().isActive()) entityManager.getTransaction().rollback();
			entityManager.close();
		} catch (final Exception exception) {
			Logger.getGlobal().log(Level.WARNING, "Closing of entity manager failed: application={0}, session={1}, message={2}.", new Object[] { session.getServletContext().getContextPath(), session.getId(), exception.getMessage() });
		}
	}


	/**
	 * Starts a JPA request demarcated "physical" transaction.
	 * @param event the request event
	 */
	public void requestInitialized (final ServletRequestEvent event) {
		if (event.getServletRequest() instanceof HttpServletRequest) {
			final HttpSession session = ((HttpServletRequest) event.getServletRequest()).getSession();
			final EntityManager entityManager = (EntityManager) session.getAttribute(JPA_ENTITY_MANAGER_KEY);
			try {
				entityManager.getTransaction().begin();
			} catch (final Exception exception) {
				Logger.getGlobal().log(Level.SEVERE, "Begin of JPA transaction failed: application={0}, session={1}, message={2}.", new Object[] { session.getServletContext().getContextPath(), session.getId(), exception.getMessage() });
				throw new IllegalStateException(exception);
			}
		}
	}


	/**
	 * Rolls back the JPA request demarcated "physical" transaction.
	 * @param event the request event
	 */
	public void requestDestroyed (final ServletRequestEvent event) {
		if (event.getServletRequest() instanceof HttpServletRequest) {
			final HttpSession session = ((HttpServletRequest) event.getServletRequest()).getSession();
			final EntityManager entityManager = (EntityManager) session.getAttribute(JPA_ENTITY_MANAGER_KEY);
			try {
				if (entityManager.getTransaction().isActive()) entityManager.getTransaction().rollback();
			} catch (final Exception exception) {
				Logger.getGlobal().log(Level.WARNING, "Rollback of JPA transaction failed: application={0}, session={1}, message={2}.", new Object[] { session.getServletContext().getContextPath(), session.getId(), exception.getMessage() });
			}
		}
	}
}