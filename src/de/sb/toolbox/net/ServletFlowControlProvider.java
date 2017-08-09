package de.sb.toolbox.net;

import java.io.IOException;
import java.net.URI;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.servlet.Filter;
import javax.servlet.FilterChain;
import javax.servlet.FilterConfig;
import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import de.sb.toolbox.Copyright;


/**
 * Generic filter class for flow control based on the phase processing design pattern. The pattern allows flexible handling of
 * page branching and input processing. Note that this flow control pattern requires that rendering is handled exclusively by a
 * separate rendering servlet or JSP. <br />
 * This implementation uses the "referer" request header field passed with POST requests to initiate a 3-phase flow:
 * <ul>
 * <li>{@code post-processing}, which processes input data and is handled by a controller</li>
 * <li>{@code pre-processing}, which processes output data and is handled by a controller</li>
 * <li>{@code rendering}, which is handled by a JSP.</li>
 * </ul>
 * The advantage is that all submit buttons and links requiring post-processing can point forward to the (projected) destination
 * page, which in turn ensures consistent web browser book-marks. However, the disadvantage is that this design requires some
 * form of controller configuration, and that an HTTP redirect is necessary whenever the server decides for a destination page
 * differing from the HTTP target.<br />
 * This filter is designed to filter all HTTP requests that it is configured for within a webapp deployment descriptor, and
 * whose path matches one of the filter's request controller dispatch paths. Requests to other paths are ignored. If a request
 * is filtered, the filter directs it's request controllers to perform phase processing. If one of the request controllers
 * involved overrules the original request target, an HTTP redirect is performed and all further actions registered down in the
 * filter chain are skipped.<br />
 * The filter may be configured within the webapp deployment descriptor to intercept all kinds of requests using the prefix
 * mapping "/*". Alternatively, it can be configured for extension mapping to intercept only requests for certain view JSPs
 * using "*.jsp" and/or "*.jspx". The following filter initialization parameters can be specified within the deployment
 * descriptor:
 * <ul>
 * <li>{@code REQUEST_ENCODING}: The encoding of incoming strings, default is UTF-8.</li>
 * <li>{@code REQUEST_CONTROLLERS}: A comma separated list of fully qualified request controller class names. The class names
 * must resolve to default constructible implementors of the PhaseProcessorServlet.RequestController interface. The first class
 * name is treated as the application's fall-back controller. Default is an empty list, but this causes all requests to be
 * bypassed by the filter.</li>
 * </ul>
 * Note that this implementation is based on servlet specification v2.5, SRV 6.2.1, stating that a container will only
 * instantiate a single instance of each configured filter per webapp context. Therefore, filter instance variables are suitable
 * for configuration information.
 */
@Copyright(year = 2005, holders = "Sascha Baumeister")
public final class ServletFlowControlProvider implements Filter {
	static private final String REQUEST_ENCODING_KEY = "REQUEST_ENCODING";
	static private final String REQUEST_CONTROLLERS_KEY = "REQUEST_CONTROLLERS";

	private final Map<String,RequestController> requestControllers = new HashMap<String,RequestController>();
	private String requestEncoding = null;



	/**
	 * Request controller interface for phase processing of HTTP requests. All request controllers registered with a phase
	 * processor filter must implement this interface. Note that all request controller classes are expected to implement a
	 * static public method {@code singleton()} returning the class's sole instance, or a public default constructor!
	 */
	static public interface RequestController {

		/**
		 * Processes the request's input data, in order to modify (persistent) data or change the web page. Returns {@code this}
		 * or an alternate controller if the operation sees a need to branch to a specific page, which may force an HTTP
		 * redirect to notify the web-client of a destination change. Alternatively returns {@code null} to indicate that the
		 * operation sees no need for branching. Note that neither case causes another input processing phase.
		 * @param request the request
		 * @return a request controller, or {@code null}
		 */
		RequestController processInputData (HttpServletRequest request);


		/**
		 * Processes the request's output data, in order to prepare for the following render phase. Returns {@code this} if
		 * processing shall continue with the rendering phase. Alternatively, returns an alternative controller if
		 * pre-processing sees a need to branch to a specific page, which may force an HTTP redirect to notify the web-client of
		 * a destination change. Note that the latter case will cause another output data processing phase for the new
		 * destination. Note that {@code null} may be returned as a synonym for {@code this}.
		 * @param request the request
		 * @return a request controller, or {@code null}
		 */
		RequestController processOutputData (HttpServletRequest request);


		/**
		 * Returns the associated (usually JSP based) renderer path.
		 * @return an absolute renderer path
		 */
		String getRendererPath ();
	}


	/**
	 * Initializes the request controllers and the request encoding from configuration data.
	 * @param filterConfiguration the filter configuration
	 */
	public void init (final FilterConfig filterConfiguration) {
		final ServletContext application = filterConfiguration.getServletContext();

		final String requestEncoding = filterConfiguration.getInitParameter(REQUEST_ENCODING_KEY);
		this.requestEncoding = requestEncoding == null || requestEncoding.isEmpty() ? "UTF-8" : requestEncoding;

		final String requestControllersText = filterConfiguration.getInitParameter(REQUEST_CONTROLLERS_KEY);
		if (requestControllersText != null) {
			final String[] classNames = requestControllersText.split(",");
			for (int index = 0; index < classNames.length; ++index) {
				final String className = classNames[index].trim();
				try {
					final Class<?> clazz = Class.forName(className, true, Thread.currentThread().getContextClassLoader());

					RequestController requestController;
					try {
						requestController = (RequestController) clazz.getDeclaredMethod("singleton").invoke(null);
					} catch (final Exception exception) {
						requestController = (RequestController) clazz.newInstance();
					}

					final String rendererPath = requestController.getRendererPath();
					if (rendererPath == null || rendererPath.isEmpty() || rendererPath.charAt(0) != '/') {
						Logger.getGlobal().log(Level.WARNING, "Ignoring illegal request controller: application={0}, class={1}, reason={2}.", new Object[] { application.getContextPath(), className, "renderer path must start with '/'" });
					} else if (this.requestControllers.containsKey(rendererPath)) {
						Logger.getGlobal().log(Level.WARNING, "Ignoring illegal request controller: application={0}, class={1}, reason={2}.", new Object[] { application.getContextPath(), className, "renderer path must be unique" });
					} else {
						this.requestControllers.put(requestController.getRendererPath(), requestController);
						if (index == 0) this.requestControllers.put("/", requestController);
					}
					Logger.getGlobal().log(Level.CONFIG, "Configured request controller: application={0}, class={1}, rendererPath={2}.", new Object[] { application.getContextPath(), className, rendererPath });
				} catch (final Exception exception) {
					Logger.getGlobal().log(Level.WARNING, "Ignoring illegal request controller: application={0}, class={1}, reason={2}.", new Object[] { application.getContextPath(), className, exception.getMessage() });
				}
			}
		}
	}


	/**
	 * Clears the configuration information stored.
	 */
	public void destroy () {
		this.requestControllers.clear();
		this.requestEncoding = null;
	}


	/**
	 * Filters the request as part of the given filter chain. The filter sets the request character encoding to the filter's
	 * encoding. Filter processing continues by calling {@link #doFilter(HttpServletRequest, HttpServletResponse, FilterChain)}
	 * if both the request and response objects are HTTP related. Otherwise, filter processing continues within the filter
	 * chain.
	 * @param request the request
	 * @param response the response
	 * @param filterChain the filter chain
	 * @throws ServletException if a servlet related problem occurs
	 * @throws IOException if the request encoding is not supported or another I/O related problem occurs
	 */
	public void doFilter (final ServletRequest request, final ServletResponse response, final FilterChain filterChain) throws ServletException, IOException {
		request.setCharacterEncoding(this.requestEncoding);

		if (request instanceof HttpServletRequest && response instanceof HttpServletResponse) {
			this.doFilter((HttpServletRequest) request, (HttpServletResponse) response, filterChain);
		} else {
			filterChain.doFilter(request, response);
		}
	}


	/**
	 * Filters the request as part of the given filter chain. The filter performs phase processing by calling
	 * {@link #doFilter(HttpServletRequest)}. Depending on the result, it either continues filtering within the filter chain, or
	 * performs an HTTP redirect.
	 * @param request the HTTP request
	 * @param response the HTTP response
	 * @param filterChain the filter chain
	 * @throws ServletException if a servlet related problem occurs
	 * @throws IOException if an I/O related problem occurs
	 */
	private void doFilter (final HttpServletRequest request, final HttpServletResponse response, final FilterChain filterChain) throws ServletException, IOException {
		final Logger logger = Logger.getGlobal();
		final RequestController redirectController = this.doFilter(request);
		if (redirectController == null) {
			if (logger.isLoggable(Level.FINER)) logger.log(Level.FINER, "Renderer for session {0} is {1}.", new Object[] { request.getSession().getId(), request.getServletPath() });
			filterChain.doFilter(request, response);
		} else {
			final String redirectPath = redirectController.getRendererPath();
			if (logger.isLoggable(Level.FINER)) logger.log(Level.FINER, "Redirecting request for session {0} to {1}.", new Object[] { request.getSession().getId(), redirectPath });
			response.sendRedirect(redirectPath.substring(1));
		}
	}


	/**
	 * Phase processes the HTTP request as part of the given filter chain. The filter determines a pre-process controller by
	 * analyzing the request URI. Processing is skipped if no pre-process controller can be determined. The filter then
	 * determines a post-process controller by analyzing the referer URI, and post-processes the current page's input.
	 * Post-Processing is skipped if no post-process controller can be determined. Afterwards, pre-processing is performed with
	 * the pre-process controller. A redirect controller is returned if any of the two controllers requests a continuation
	 * different from the one defined by the request URI, otherwise {@code null}.
	 * @param request the HTTP request
	 * @return a redirect controller or {@code null}
	 */
	private RequestController doFilter (final HttpServletRequest request) {
		final Logger logger = Logger.getGlobal();
		final String localRendererPath = request.getServletPath(), localRefererPath = localRefererPath(request);

		RequestController nextController = this.requestControllers.get(localRefererPath);
		if (nextController != null) {
			if (logger.isLoggable(Level.FINER)) logger.log(Level.FINER, "Post-processor for session {0} is {1}.", new Object[] { request.getSession().getId(), nextController.getClass().getName() });
			nextController = nextController.processInputData(request);
			if (nextController != null && !nextController.getRendererPath().equals(localRendererPath)) return nextController;
		}

		if (nextController == null) nextController = this.requestControllers.get(localRendererPath);
		if (nextController != null) {
			if (logger.isLoggable(Level.FINER)) logger.log(Level.FINER, "Pre-processor for session {0} is {1}.", new Object[] { request.getSession().getId(), nextController.getClass().getName() });
			nextController = nextController.processOutputData(request);
			if (nextController != null && !nextController.getRendererPath().equals(localRendererPath)) return nextController;
		}

		return null;
	}


	/**
	 * Returns the local (i.e. relativized) referer path if the request's referer URL is present and correctly formatted,
	 * targets the local HTTP server, and the URL's path starts with the local application's context path. Returns {@code null}
	 * otherwise.
	 * @param request the HTTP request
	 * @return a local referer path, or {@code null}
	 */
	static private String localRefererPath (final HttpServletRequest request) {
		final URI refererURI;
		try {
			refererURI = new URI(request.getHeader("referer"));
		} catch (final Exception exception) {
			return null;
		}

		if (refererURI.getHost().equals(request.getServerName()) && refererURI.getPort() == request.getServerPort()) {
			final String contextPath = request.getContextPath();
			if (refererURI.getPath().startsWith(contextPath)) {
				return refererURI.getPath().substring(contextPath.length());
			}
		}
		return null;
	}
}