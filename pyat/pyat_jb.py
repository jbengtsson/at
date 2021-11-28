import inspect
import os, sys
import at
import at.plot
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (9.0, 6.0)


def prt_ps(str, ps):
    print(str, end='')
    for x in ps:
        print('{:11.3e}'.format(x), end='')
    print()


def prt_mat(str, a):
    print(str, end='')
    for j in range(0, len(a)):
        for k in range(0, len(a[0])):
            print('{:14.6e}'.format(a[j][k]), end='')
        print()
    

def get_optics(ring, plot):
    [elemdata0, beamdata, elemdata] = \
        at.get_optics(ring, range(len(ring)+1), get_chrom=True)
    print('\n  nu  = [{:5.3f}, {:5.3f}]'
          .format(beamdata.tune[0], beamdata.tune[1]))
    print('  ksi = [{:5.3f}, {:5.3f}]'
          .format(beamdata.chromaticity[0], beamdata.chromaticity[1]))
    [M, _] = at.find_m66(ring, 0)
    prt_mat('\nPoincaré map:\n', M)
    print('\n', at.radiation_parameters(ring))
    if plot:
        ring.plot_beta()
        plt.show()


def example(ring):
    if False:
        print('\n', dir(ring))

    if False:
        for elem in ring:
            print(elem)

    if False:
        refq1 = at.get_cells(ring, at.checktype(at.Quadrupole))
        print('\n', refq1)
        print('\n', list(ring[refq1]))

    if not False:
        ps = np.array([1e-3, -1e-3, 0e0, 0e0, 0e0, 0e0])
        prt_ps('\nTrack one turn:\n', ps)
        at.lattice_pass(sup_per, ps, 1)
        prt_ps('', ps)

    if False:
        [cod, _] = at.find_orbit(sup_per, dp=0e0)
        prt_ps('\nFixed point:\n', cod)

    if False:
        [elemdata0, beamdata, elemdata] = \
            at.get_optics(sup_per, range(len(sup_per)-1))
        print('\nelemdata.shape:', elemdata.shape)
        print('elemdata.fields:')
        for fld in elemdata.dtype.fields.keys():
            print(fld)

    if False:
        plt.plot(elemdata.s_pos, elemdata.beta)
        plt.xlabel('s [m]')
        plt.ylabel(r'$\beta$ [m]')
        plt.show()

    if not False:
        ring.plot_beta()
        plt.show()

    if not False:
        # Turn on RF cavity.
        # Assumes that it's at end of lattice.
        n = len(ring)
        print('\n', ring[n-1])
        ring[n-1].PassMethod = "CavityPass"
        print('\n', ring[n-1])

    if False:
        ring.radiation_on(ring)
        [_, beamdata, _] = at.ohmi_envelope(ring)
        print('\nbeamdata.fields:')
        for fld in beamdata.dtype.fields.keys():
            print(fld)

def tst_cases:
    import numpy
    from at import elements, atpass

    # test_missing_pass_method_raises_attribute_error.
    lat = [elements.Marker('marker')]
    ps = numpy.zeros(6,)
    print('\n', lat[0])
    del lat[0].PassMethod
    atpass(lat, ps, 1)

    # test_missing_length_raises_attribute_error.
    lat = [elements.Drift('drift', 1.0)]
    ps = numpy.zeros(6,)
    ps[1] = 1e-6
    del lat[0].Length
    atpass(lat, ps, 1)

    # test_quad.
    q = elements.Quadrupole('quad', 0.4, k=1)
    lattice = [q]
    ps = numpy.array(numpy.zeros((6, 1)), order='F')
    ps[0, 0] = 1e-6
    atpass(lattice, ps, 1)
    print(ps)

    # test_rfcavity.
    rf = elements.RFCavity('rfcavity', 0.0, 187500, 3.5237e+8, 31, 6.e+9)
    print(rf)
    lattice = [rf, rf, rf, rf]
    ps = numpy.array(numpy.zeros((6, 1)), order='F')
    ps[4, 0] = 1e-6
    ps[5, 0] = 1e-6
    atpass(lattice, ps, 1)
    print(ps)

lat_dir = os.environ['LAT']

lat_name = {
    'hmba'      : 'test_matlab/hmba.mat',
    'esrf_ebs'  : 'test_matlab/err.mat',
    'bessy-iii' : lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy.lat'
}

lat = 'bessy-iii'
if lat != 'bessy-iii':
    sup_per = at.load_mat(lat_name[lat])
else:
    sup_per = \
        at.load_tracy(lat_name[lat], harmonic_number=538, lattice_key='cell')

if False:
    for elem in sup_per:
        print(elem)

#print('\n', sup_per)

sup_per.radiation_off()

n = len(sup_per)
if False:
    # Turn on RF cavity.
    # Assumes that it's at end of lattice.
    print('\n', sup_per[n-1])
    sup_per[n-1].PassMethod = "CavityPass"
    print('\n', sup_per[n-1])

#example(sup_per)

get_optics(sup_per, False)


# AT Interface
#
# Files:
#   ../at/atintegrators/*.m *.c *.h
#   ../at/pyat/at.c
#             /at/integrators/*.so
#             /integrator-src/*.h *.c
#
# Propagators:
#
#   IdentityPass()
#   DriftPass()
#   StrMPoleSymplectic4Pass()
#   BndMPoleSymplectic4Pass()
#   CavityPass()

# gen_m66_elem(), gen_detuning_elem(), gen_quantdiff_elem()
# element_pass(), lattice_pass()
# find_orbit(), find_orbit4(), find_orbit6()
# find_m44(), find_m66()
# get_tune()
# linopt(), linopt6()
#
# [elemdata0, beamdata, elemdata] = \
#     at.get_optics(ring, range(len(ring)+1), get_chrom=True)
#
# [M, _] = at.find_m66(ring, 0)
#
# print('\n', at.radiation_parameters(ring))
#
# at.lattice_pass(sup_per, ps, 1)
#
# [_, beamdata, _] = at.ohmi_envelope(ring)
