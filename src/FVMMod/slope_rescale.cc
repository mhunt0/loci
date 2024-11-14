/** ****************************************************************************
 * @file      slope_rescale.cc
 * @author    Ed Luke
 * @date      LICENSE Date: 12-30-2023
 * @copyright MS State/CFDRC
 * @brief     Sets frozen/minimum limiter for vect3d state variables
 * @details   The software tool Loci/GGFS module consisting is being furnished
 *            to NASA Personnel under SBIR Data Rights. The SBIR DATA rights are
 *            asserted by CFD Research Corporation.
 *
 *            These SBIR data are furnished with SBIR rights under Contract No.
 *            80NSSC18P2154. For a period of 4 years, unless extended in
 *            accordance with FAR 27.409(h), after acceptance of all items to be
 *            delivered under this contract, the Government will use these data
 *            for Government purposes only, and they shall not be disclosed
 *            outside the Government (including disclosure for procurement
 *            purposes) during such period without permission of the Contractor,
 *            except that, subject to the foregoing use and disclosure
 *            prohibitions, these data may be disclosed for use by support
 *            Contractors. After the protection period, the Government has a
 *            paid-up license to use, and to authorize others to use on its
 *            behalf, these data for Government purposes, but is relieved of all
 *            disclosure prohibitions and assumes no liability for unauthorized
 *            use of these data by third parties. This notice shall be affixed
 *            to any reproductions of these data, in whole or in part.
 * @attention Distribution C: Limited to Government Employees only. A release
 *            under Distribution B and A is being considered and may be done for
 *            future releases of the code.
 ******************************************************************************/

#include "slope_rescale.h"
    
namespace Loci {
/** ****************************************************************************
 * @brief Void function to rescale slope in muscl reconstruction for vectors
 * @param[out] slopes  real, rescaled slope for reconstuction
 * @param[in]  vs      int, size of vector
 ******************************************************************************/
    void slope_rescale(real *slopes, int vs) {

      real sum = 0.0, sumn = 0.0, sump = 0.0 ;
      for(int i=0;i<vs;++i) {
	sum += slopes[i] ;
	sumn += min(real(0.0),slopes[i]) ;
	sump += max(real(0.0),slopes[i]) ;
      }
      if(sum < 0.0) {
	for(int i=0;i<vs;++i)
	  slopes[i] = slopes[i]*((slopes[i]<0.0)?(sumn-sum)/(sumn-1e-30):1.0) ;
      }
      if(sum > 0.0) {
	for(int i=0;i<vs;++i)
	  slopes[i] = slopes[i]*((slopes[i]>0.0)?(sump-sum)/(sump-1e-30):1.0) ;
      }
    }    
}