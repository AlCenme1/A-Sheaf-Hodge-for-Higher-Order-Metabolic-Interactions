"""Stage 1-2: SBML parsing and organism/metabolite identification."""
import xml.etree.ElementTree as ET
import re, json, time

NS = {'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
      'fbc': 'http://www.sbml.org/sbml/level3/version1/fbc/version2'}

def parse_sbml(path):
    t0=time.time()
    tree = ET.parse(path)
    root = tree.getroot()
    model = root.find('sbml:model', NS)

    species = {}
    for sp in model.find('sbml:listOfSpecies', NS):
        species[sp.get('id')] = {'compartment': sp.get('compartment'), 'name': sp.get('name',''),
                                  'boundary': sp.get('boundaryCondition', 'false') == 'true'}

    reactions = {}
    for rxn in model.find('sbml:listOfReactions', NS):
        rid = rxn.get('id')
        lb_ref = rxn.get('{%s}lowerFluxBound' % NS['fbc'])
        ub_ref = rxn.get('{%s}upperFluxBound' % NS['fbc'])
        reactants, products = {}, {}
        lor = rxn.find('sbml:listOfReactants', NS)
        if lor is not None:
            for sr in lor:
                reactants[sr.get('species')] = float(sr.get('stoichiometry', 1))
        lop = rxn.find('sbml:listOfProducts', NS)
        if lop is not None:
            for sp_ in lop:
                products[sp_.get('species')] = float(sp_.get('stoichiometry', 1))
        reactions[rid] = {'reactants': reactants, 'products': products,
                           'lb_ref': lb_ref, 'ub_ref': ub_ref}

    params = {}
    lop_params = model.find('sbml:listOfParameters', NS)
    if lop_params is not None:
        for p in lop_params:
            params[p.get('id')] = float(p.get('value', 0))
    print(f"parsed in {time.time()-t0:.1f}s: {len(species)} species, {len(reactions)} reactions, {len(params)} params")
    return species, reactions, params

if __name__ == '__main__':
    from paths import COMMUNITY
    species, reactions, params = parse_sbml(str(COMMUNITY))
    print(len(species), 'species,', len(reactions), 'reactions')
