/*
mutate
   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "mutate.h"

sz count_mutable_entities(const conf& c) {
	sz counter = 0;
	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();
	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	return counter;
}

// does not set model
void mutate_conf(conf& c, const model& m, fl amplitude, rng& generator) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	sz mutable_entities_num = count_mutable_entities(c);
	if(mutable_entities_num == 0) return;
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	VINA_CHECK(which_int >= 0);
	sz which = sz(which_int);
	VINA_CHECK(which < mutable_entities_num);

	VINA_FOR_IN(i, c.ligands) {
		if(which == 0) { 
			std::cout << "mutate position" << std::setw(15);
			c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator);
			return; 
		}
		--which;
		if(which == 0) { 
			std::cout << "mutate orientation" << std::setw(15);
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
				vec rotation; 
				rotation = amplitude / gr * random_inside_sphere(generator); 
				quaternion_increment(c.ligands[i].rigid.orientation, rotation);
			}
			return; 
		}
		--which;
		if(which < c.ligands[i].torsions.size()) { 
			std::cout << "mutate torsion" <<  which + 1 << std::setw(15);
			c.ligands[i].torsions[which] = random_fl(-pi, pi, generator); 
			return; 
		}
		which -= c.ligands[i].torsions.size();
	}
	VINA_FOR_IN(i, c.flex) {
		if(which < c.flex[i].torsions.size()) { c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= c.flex[i].torsions.size();
	}
}
 
void mutate_conf_partial(conf& c, const model& m, fl amplitude, rng& generator) {
	fl mutate_pro = 0.5;
	fl which = random_fl(0, 1, generator);

	VINA_FOR_IN(i, c.ligands) {
		if (which <= mutate_pro) {
			std::cout << "mutate orientation" << std::setw(15);
			fl gr = m.gyration_radius(i);
			if (gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
				vec rotation;
				rotation = amplitude / gr * random_inside_sphere(generator);
				quaternion_increment(c.ligands[i].rigid.orientation, rotation);
			}
			return;
		}
		else if (which > mutate_pro) {
			std::cout << "mutate position" << std::setw(15);
			c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator); 
			return;
		}
	}
}
