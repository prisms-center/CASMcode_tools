"""Search for mappings between two structures."""

import math

import casm.map.utils as utils
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal


def search(args):
    """
    This search has the following behavior:
    - reads a parent structure and child structure
    - optionally masks occupant types, treating all species as equal
    - optionally allows vacancies on the parent structure
    - optionally makes the parent or the child or both primitive
    - optionally constructs supercells of the child up to a certain volume
    - constructs the factor groups of the parent and child
    - checks structure compatibility as detailed below
    - performs a mapping with max_vol, k_best, min_cost, max_cost, cost type
      - defaults: max_vol determined automatically, rest set to map_structures defaults
      - max_vol default is lcm(child's # atoms, parent's # sites) / (parent's # sites)
      - if n vacancies are allowed, it is lcm(child's # atoms + n, parent's # sites) /
        (parent's # sites)
    - returns mapping results and the parent and child structures which were used
    - note: vacancy fraction needs to be handled somehow, in case the user passes in
      weird numbers

    This search performs the following checks:
    - the child's number of atoms must be >= the parent's number of sites
      - if vacancies are allowed, (child's # atoms + # Va allowed) >= (parent's # sites)
    - the child's number of atoms must be an integer multiple of the parent's number of
      sites
      - if max_vol is unset (default), the max_vol is set automatically and
        lcm(child's # atoms, parent's # sites) == (child's # atoms) must be true
      - if vacancies are allowed, the following must be true for some n in
        [0, # Va allowed]:
        lcm(child's # atoms + n, parents # sites) == (child's # atoms + n)
    - the occupants must be reasonable (TODO)
    - k_best must be >=1
    """
    if args.k_best < 1:
        raise ValueError(f"k_best must be >=1 (got {args.k_best})")
    k_best = args.k_best

    if args.child_supercells < 1:
        raise ValueError(f"child_supercells must be >= 1 (got {args.child_supercells})")

    # read parent and child
    primify_parent = "parent" in args.primify or "both" in args.primify
    # symmetrize_parent = "parent" in args.symmetrize or "both" in args.symmetrize
    primify_child = "child" in args.primify or "both" in args.primify
    # symmetrize_child = "child" in args.symmetrize or "both" in args.symmetrize
    print(f"reading parent: {args.parent}, primify: {primify_parent}")
    parent = utils.read_prim(args.parent, primify=primify_parent, symmetrize=False)
    print(f"reading child: {args.child}, primify: {primify_child}")
    child = utils.read_structure(args.child, primify=primify_child, symmetrize=False)

    # should come before include_vacancies
    if args.mask_occupants:
        parent_dict = parent.to_dict()
        parent_dict["basis"] = [
            {"coordinate": site["coordinate"], "occupants": ["X"]}
            for site in parent_dict["basis"]
        ]
        parent = xtal.Prim.from_dict(parent_dict)
        child_dict = child.to_dict()
        child_dict["atom_type"] = ["X"] * len(child_dict["atom_type"])
        child = xtal.Structure.from_dict(child_dict)

    # get max volume for parent
    # if args.include_vacancies:
    #     raise NotImplementedError
    # if args.child_supercells:
    #     raise NotImplementedError
    if args.max_vol is None:
        parent_natoms = parent.coordinate_cart().shape[1]
        child_natoms = len(child.atom_type())
        if parent_natoms > child_natoms:
            print(
                f"incompatible structures: parent has more atoms ({parent_natoms}) "
                f"than child ({child_natoms})\n"
                "try making supercells of the child or primifying the parent"
            )
            return
        if math.lcm(parent_natoms, child_natoms) != child_natoms:
            print(
                f"incompatible structures: child natoms ({child_natoms}) is not an "
                f"integer multiple of parent natoms ({parent_natoms})\n"
                f"try making supercells of the structure or primifying them"
            )
            return
        max_vol = math.lcm(parent_natoms, child_natoms) // parent_natoms
    else:
        max_vol = args.max_vol

    # make factor groups
    # if len(args.symmetrize_if_necessary) != 0:
    #     raise NotImplementedError
    prim_factor_group = xtal.make_factor_group(parent)
    structure_factor_group = xtal.make_factor_group(child)

    # map
    ## note: will have to be changed if multiple supercells are created
    maps = mapmethods.map_structures(
        prim=parent,
        structure=child,
        max_vol=max_vol,
        prim_factor_group=prim_factor_group,
        structure_factor_group=structure_factor_group,
        k_best=k_best,
    )
    return maps, parent, child


def main():
    print("search")


if __name__ == "__main__":
    main()


def search_tmp(child, args):
    """
    temp search method that takes a child instead of a path to one
    """
    if args.k_best < 1:
        raise ValueError(f"k_best must be >=1 (got {args.k_best})")
    k_best = args.k_best

    if args.child_supercells < 1:
        raise ValueError(f"child_supercells must be >= 1 (got {args.child_supercells})")

    # read parent and child
    primify_parent = "parent" in args.primify or "both" in args.primify
    # symmetrize_parent = "parent" in args.symmetrize or "both" in args.symmetrize
    #    primify_child = "child" in args.primify or "both" in args.primify
    # symmetrize_child = "child" in args.symmetrize or "both" in args.symmetrize
    print(f"reading parent: {args.parent}, primify: {primify_parent}")
    parent = utils.read_prim(args.parent, primify=primify_parent, symmetrize=False)

    # should come before include_vacancies
    if args.mask_occupants:
        parent_dict = parent.to_dict()
        parent_dict["basis"] = [
            {"coordinate": site["coordinate"], "occupants": ["X"]}
            for site in parent_dict["basis"]
        ]
        parent = xtal.Prim.from_dict(parent_dict)
        child_dict = child.to_dict()
        child_dict["atom_type"] = ["X"] * len(child_dict["atom_type"])
        child = xtal.Structure.from_dict(child_dict)

    # get max volume for parent
    # if args.include_vacancies:
    #     raise NotImplementedError
    # if args.child_supercells:
    #     raise NotImplementedError
    if args.max_vol is None:
        parent_natoms = parent.coordinate_cart().shape[1]
        child_natoms = len(child.atom_type())
        if parent_natoms > child_natoms:
            print(
                f"incompatible structures: parent has more atoms ({parent_natoms}) "
                f"than child ({child_natoms})\n"
                "try making supercells of the child or primifying the parent"
            )
            return
        if math.lcm(parent_natoms, child_natoms) != child_natoms:
            print(
                f"incompatible structures: child natoms ({child_natoms}) is not an "
                f"integer multiple of parent natoms ({parent_natoms})\n"
                f"try making supercells of the structure or primifying them"
            )
            return
        max_vol = math.lcm(parent_natoms, child_natoms) // parent_natoms
    else:
        max_vol = args.max_vol

    # make factor groups
    # if len(args.symmetrize_if_necessary) != 0:
    #     raise NotImplementedError
    prim_factor_group = xtal.make_factor_group(parent)
    structure_factor_group = xtal.make_factor_group(child)

    # map
    ## note: will have to be changed if multiple supercells are created
    maps = mapmethods.map_structures(
        prim=parent,
        structure=child,
        max_vol=max_vol,
        prim_factor_group=prim_factor_group,
        structure_factor_group=structure_factor_group,
        k_best=k_best,
    )
    return maps, parent, child
