package de.sb.toolbox.net;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import javax.persistence.EntityManager;
import javax.persistence.EntityManagerFactory;
import javax.persistence.Persistence;
import javax.ws.rs.container.ContainerRequestContext;
import javax.ws.rs.container.ContainerRequestFilter;
import javax.ws.rs.container.ContainerResponseContext;
import javax.ws.rs.container.ContainerResponseFilter;
import javax.ws.rs.ext.Provider;
import de.sb.toolbox.Copyright;


/**
 * This lifecycle provider deploys the following services for a given JPA persistence unit:
 * <ul>
 * <li><b>Entity manager life-cycle management</b>: Entity managers for this provider's persistence unit are created upon the
 * begin of any HTTP request, and closed upon it's end. This implies that the entity managers are active during entity
 * marshaling, ready to supply additional information.</li>
 * <li><b>Request-scoped transaction demarcation</b>: Additionally, the design allows for continuous transaction coverage,
 * similar to JDBC. The idea is that a transaction is started automatically upon request begin, and at it's end the last active
 * transaction is automatically rolled back. Services should immediately start a new transaction after committing an existing
 * one.</li>
 * </ul>
 * Note that the use of a thread local variable for entity manager injection is based on the precondition that any HTTP request
 * is processed within a single thread. This assumption does hold in standard compatible environments, like Jersey.
 */
@Provider
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class RestJpaLifecycleProvider implements ContainerRequestFilter, ContainerResponseFilter, AutoCloseable {
	static private final Map<String,ThreadLocal<EntityManager>> THREAD_LOCAL_ENTITY_MANAGERS = Collections.synchronizedMap(new HashMap<>());

	private final String persistenceUnitName;
	private final EntityManagerFactory entityManagerFactory;
	private final ThreadLocal<EntityManager> entityManagerReference;


	/**
	 * Returns the entity manager associated with both the current thread and the given persistence unit.
	 * @param persistenceUnitName the persistence unit name
	 * @return the entity manager
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if there is no lifecycle provider associated with the given persistence unit
	 * @throws IllegalStateException if there is no entity manager associated with the current thread
	 */
	static public EntityManager entityManager (final String persistenceUnitName) throws NullPointerException, IllegalArgumentException, IllegalStateException {
		if (persistenceUnitName == null) throw new NullPointerException();

		final ThreadLocal<EntityManager> entityManagerReference = THREAD_LOCAL_ENTITY_MANAGERS.get(persistenceUnitName);
		if (entityManagerReference == null) throw new IllegalArgumentException();

		final EntityManager entityManager = entityManagerReference.get();
		if (entityManager == null) throw new IllegalStateException();

		return entityManager;
	}


	/**
	 * Creates a new instance, and forces entity manager factory initialization if it has not happened yet. Note that this
	 * allows JAX-RS to create multiple instances sharing the same entity manager factory.
	 * @param persistenceUnitName the persistence unit name
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument differs from the name passed with a previous instantiation
	 * @throws RuntimeException if there is a problem configuring the persistence unit
	 */
	public RestJpaLifecycleProvider (final String persistenceUnitName) throws RuntimeException {
		if (persistenceUnitName == null) throw new NullPointerException();

		this.persistenceUnitName = persistenceUnitName;
		this.entityManagerFactory = Persistence.createEntityManagerFactory(persistenceUnitName);
		this.entityManagerReference = new ThreadLocal<>();

		synchronized (THREAD_LOCAL_ENTITY_MANAGERS) {
			if (THREAD_LOCAL_ENTITY_MANAGERS.containsKey(this.persistenceUnitName)) throw new IllegalArgumentException();
			THREAD_LOCAL_ENTITY_MANAGERS.put(this.persistenceUnitName, this.entityManagerReference);
		}
	}


	/**
	 * {@inheritDoc}
	 */
	public void close () {
		try {
			THREAD_LOCAL_ENTITY_MANAGERS.remove(this.persistenceUnitName);
		} finally {
			this.entityManagerFactory.close();
		}
	}


	/**
	 * Ensures that this provider is closed before the garbage collector removes it. This is not strictly neccessary, but may
	 * free resources in a more timely manner.
	 */
	@Override
	protected void finalize () {
		this.close();
	}


	/**
	 * This operation is called by the JAX-RS runtime before an HTTP request is processed withing the current thread. It creates
	 * a new entity manager instance using this provider's entity manager factory, and stores it within the thread local
	 * reference associated with this provider's persistence unit name.
	 * @param requestContext the (optional) JAX-RS request context
	 */
	public void filter (final ContainerRequestContext requestContext) {
		final EntityManager entityManager = this.entityManagerFactory.createEntityManager();
		entityManager.getTransaction().begin();
		this.entityManagerReference.set(entityManager);
	}


	/**
	 * This operation is called by the JAX-RS runtime after an HTTP request has been processed, but before the entity stream has
	 * been written. It decorates the response context's entity stream, causing said decorator to trigger once the entity stream
	 * has been written. This in turn allows this operation to close and remove the entity manager associated with both the
	 * current thread and this provider's persistence unit. Note that this technology relies on the entity stream (rather, the
	 * decorator wrapping it) to be closed regardless of the presence of absence of a response entity; in other words, the
	 * operation relies heavily on correct resource management by the JAX-RS implementation.
	 * @param requestContext the JAX-RS request context
	 * @param responseContext the JAX-RS response context
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 */
	public void filter (final ContainerRequestContext requestContext, final ContainerResponseContext responseContext) throws NullPointerException {
		final FilterOutputStream triggerStream = new FilterOutputStream(responseContext.getEntityStream()) {

			/**
			 * {@inheritDoc}
			 */
			public void close () throws IOException {
				try {
					super.close();
				} finally {
					RestJpaLifecycleProvider.this.closeEntityManager();
				}
			}
		};
		responseContext.setEntityStream(triggerStream);
	}


	/**
	 * Closes the entity manager associated with both the current thread and this provider's persistence unit.
	 */
	private void closeEntityManager () {
		final EntityManager entityManager = this.entityManagerReference.get();
		this.entityManagerReference.remove();

		if (entityManager != null && entityManager.isOpen()) {
			try {
				if (entityManager.getTransaction().isActive()) entityManager.getTransaction().rollback();
			} finally {
				entityManager.close();
			}
		}
	}
}