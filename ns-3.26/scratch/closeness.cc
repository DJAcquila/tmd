/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "ns3/core-module.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

using namespace ns3;


double my_integrand (double x, void *params);


NS_LOG_COMPONENT_DEFINE ("ScratchSimulator");

int 
main (int argc, char *argv[])
{
  NS_LOG_UNCOND ("Scratch Simulator");

 gsl_integration_workspace *work_ptr
    = gsl_integration_workspace_alloc (1000);

  double lower_limit = 0;	/* lower limit a */
  double upper_limit = 20.0/4.0;/* upper limit b */
  double abs_error = 1.0e-8;	/* to avoid round-off problems */
  double rel_error = 1.0e-8;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  double alpha = 7.0;		// parameter in integrand
  //double expected = -4.0;	// exact answer

  gsl_function My_function;
  void *params_ptr = &alpha;

  My_function.function = &my_integrand;
  My_function.params = params_ptr;

  gsl_integration_qags (&My_function, lower_limit, upper_limit,
			abs_error, rel_error, 1000, work_ptr, &result, &error);

  std::cout.precision (18);  
  //int width = 20;  // setw width for output
  std::cout << "result          = " << result << std::endl;
  std::cout << "estimated error = " << error << std::endl;

  double gamma = gsl_sf_gamma (alpha);

  std::cout << "gamma          = " << gamma << std::endl;

  double closeness = 1 - (result/gamma);

  std::cout << "closeness (wij) = " << closeness << std::endl;


  Simulator::Run ();
  Simulator::Destroy ();
}

double
my_integrand (double x, void *params)
{
  // Mathematica form: (t^(k-1))*(exp(-t))

  // The next line recovers alpha from the passed params pointer
  double alpha = *(double *) params;

  return ((pow(x,(alpha -1))) * (exp(-x)) );

}
