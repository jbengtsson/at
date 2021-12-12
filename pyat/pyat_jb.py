import inspect
import os, sys
import at
import at.plot
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (9.0, 6.0)


def printf(format, *args): sys.stdout.write(format % args)


def prt_ps(str, ps):
    print(str, end='')
    for x in ps:
        print('{:11.3e}'.format(x), end='')
    print()


def prt_mat(str, a):
    print(str, end='')
    for j in range(0, len(a)):
        for k in range(0, len(a[0])):
            if True:
                print('{:14.6e}'.format(a[j][k]), end='')
            else:
                print('{:23.15e}'.format(a[j][k]), end='')
        print()
    

def prt_elem(str, elem, verbose):
    if verbose:
        print(str, elem)
    else:
        printf('%s %-10s %10s', str, type(elem).__name__, elem.FamName)


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
        at.lattice_pass(lat, ps, 1)
        prt_ps('', ps)

    if False:
        [cod, _] = at.find_orbit(lat, dp=0e0)
        prt_ps('\nFixed point:\n', cod)

    if False:
        [elemdata0, beamdata, elemdata] = \
            at.get_optics(lat, range(len(lat)-1))
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

def tst_cases(lat):
    import numpy
    from at import elements, atpass

    if False:
        # test_missing_pass_method_raises_attribute_error.
        lat = [elements.Marker('marker')]
        ps = numpy.zeros(6,)
        print('\n', lat[0])
        del lat[0].PassMethod
        # Lattice, phase-space vector, numer of particles.
        atpass(lat, ps, 1)

    if False:
        # test_missing_length_raises_attribute_error.
        lat = [elements.Drift('drift', 1.0)]
        ps = numpy.zeros(6,)
        ps[1] = 1e-6
        del lat[0].Length
        atpass(lat, ps, 1)
        exit(0)

    if False:
        # test_quad.
        q = elements.Quadrupole('quad', 0.4, k=1)
        lattice = [q]
        ps = numpy.array(numpy.zeros((6, 1)), order='F')
        ps[0, 0] = 1e-6
        atpass(lat, ps, 1)
        print('\n', ps)
        exit(0)

    if False:
        # test_rfcavity.
        rf = elements.RFCavity('rfcavity', 0.0, 187500, 3.5237e+8, 31, 6.e+9)
        print('\n', rf)
        lattice = [rf, rf, rf, rf]
        ps = numpy.array(numpy.zeros((6, 1)), order='F')
        ps[4, 0] = 1e-6
        ps[5, 0] = 1e-6
        atpass(lattice, ps, 1)
        print('\n', ps)
        exit(0)

    if False:
        # Compare the results of linopt2 and linopt6 in 4d.
        from at import linopt2, linopt6
        from numpy.testing import assert_allclose as assert_close

        if True:
            lat.radiation_off()
        else:
            lat.radiation_on()

        print('\n', lat)

        refpts = range(len(lat) + 1)
        print('\nrefpts:\n', refpts)
        ld02, rd2, ld2 = linopt2(lat, refpts, get_w=True)
        print('\nlinopt2:\n', ld02)
        ld06, rd6, ld6 = linopt6(lat, refpts, get_w=True)
        print('\nlinopt6\n', ld06)

        # rd2.tune
        # rd2.chromaticity

        # assert_close(rd2.tune, rd6.tune, atol=1e-12, rtol=0)
        # assert_close(rd2.chromaticity, rd6.chromaticity, atol=1e-12, rtol=0)

        for field in ['s_pos', 'closed_orbit', 'dispersion', 'alpha', 'beta',
                      'mu']:
            assert_close(ld2[field], ld6[field], atol=1e-10, rtol=0,
                         err_msg=field)
            assert_close(ld2.W, ld6.W, atol=1e-6, rtol=0)
        exit(0)

    if not False:
        lat.radiation_off()

        print('\n', lat)

        ref_pts = range(len(lat)+1)
        print('\nref_pts:\n', ref_pts)

        if False:
            [M, Ms] = at.find_m66(lat, ref_pts)
            for k in range(Ms.shape[0]-1):
                printf('\n%4d', k)
                prt_elem('', lat[k], False)
                prt_mat('\n', Ms[k,:,:])

        [orbit, _] = at.find_orbit4(lat, dp=0e0)
        prt_ps('\norbit: ', orbit)

        ld02, rd2, ld2 = at.linopt2(lat, ref_pts, get_w=True)
        print('\nlinopt2:\n', ld02)
        ld06, rd6, ld6 = at.linopt6(lat, ref_pts, get_w=True)
        print('\nlinopt6\n', ld06)
        exit(0)


lat_dir = os.environ['LAT']

lat_names = {
    # 'dba'       : 'test_matlab/dba.mat',
    'dba'       : '/Users/johan/git/at/pyat/test_matlab/dba.mat',
    'hmba'      : 'test_matlab/hmba.mat',
    'esrf_ebs'  : 'test_matlab/err.mat',
    'bessy-iii' : lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy.lat'
}

# lat_name = 'dba'
lat_name = 'hmba'
# lat_name = 'bessy-iii'
if lat_name != 'bessy-iii':
    lat = at.load_mat(lat_names[lat_name])
else:
    lat = \
        at.load_tracy(lat_names[lat_name], harmonic_number=538,
                      lattice_key='cell')

if False:
    for elem in lat:
        print(elem)
    exit(0)

lat.radiation_off()

if False:
    tst_cases(lat)
    exit(0)

#print('\n', lat)

n = len(lat)
if False:
    # Turn on RF cavity.
    # Assumes that it's at end of lattice.
    print('\n', lat[n-1])
    lat[n-1].PassMethod = "CavityPass"
    print('\n', lat[n-1])

#example(lat)

get_optics(lat, False)


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
# at.lattice_pass(lat, ps, 1)
#
# [_, beamdata, _] = at.ohmi_envelope(ring)
