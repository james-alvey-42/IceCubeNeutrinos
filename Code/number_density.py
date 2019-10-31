from integration import integrate_trap
import numpy as np

def integrand(z, OmegaM0, OmegaL0):
	r"""
	Returns the integrand for the redshift of the number density,

	.. math:: \frac{(1 + z)^2}{\sqrt{\Omega_\Lambda + \Omega_{m, 0}(1 + z)^3}}

	Parameters
	----------
	z : float
		redshift
	OmegaMO : float
		matter desnity parameter today
	OmegaL0 : float
		dark energy density parameter

	Returns
	-------
	integrand : float
		integrand as given above

	Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from number_density import integrand;

       In [2]: integrand(z=0.36, OmegaM0=0.315, OmegaL0=0.65)
	"""
	numerator = np.power(1.0 + z, 2)
	denominator = np.power(OmegaL0 + OmegaM0*np.power(1.0 + z, 3), 0.5)
	return numerator * np.power(denominator, -1)

if __name__ == '__main__':
	zb = 0.3365
	c = 3.0*np.power(10.0, 5)
	OmegaM = 0.315
	OmegaL = 0.65
	H0 = 6.73*np.power(10.0, 4)
	z_arr = np.linspace(0, zb, 1000)
	int_arr = integrand(z_arr, OmegaM, OmegaL)
	integration = integrate_trap(z_arr, int_arr)
	integration = c*np.power(H0, -1)*integration
	print(integration)
