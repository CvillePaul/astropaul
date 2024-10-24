import math
import numpy as np
import pandas as pd


def get_object_flux_dictionary():

    spectypes = ["O5V", "O9V", "B0V", "B1V", "B3V", "B6V", "B8V", "A0V", "A2V", "A3V", "A5V", "F0V", "F2V"]
    spec_str = "F5V, F8V, G0V, G2V, G5V, G8V, K0V, K2V, K5V, K7V, M0V, M2V, M4V, M5V, B2IV, B6IV, A0IV, "
    spec_str = spec_str + "A4-7IV, F0-2IV, F5IV, F8IV, G0IV, G2IV, G5IV, G8IV, K0IV, K1IV, K3IV, O8III, B1-2III, B5III, "
    spec_str = spec_str + "B9III, A0III, A5III, F0III, F5III, G0III, G5III, G8III, K0III, K3III, K5III, M0III, M5III, "
    spec_str = spec_str + "M10III, B2II, B5II, F0II, F2II, G5II, K0-1II, K3-4II, M3II, B0I, B5I, B8I, A0I, F0I, "
    spec_str = spec_str + "F5I, F8I, G0I, G5I, G8I, K2I, K4I, M2I, Seyfert1, Seyfert2, QSO, NGC1068, Liner"

    spectypes = spectypes + spec_str.split(", ")

    flux_array = [
        [239593039.763, 178543473.813, 124031554.257, 82534178.9555, 52094588.3886, 30216798.7116],
        [224137694.103, 172205034.458, 126093669.89, 80045933.2903, 51523158.0048, 30126736.1341],
        [228955744.058, 174893335.89, 123225191.932, 83093894.2227, 53409615.9756, 31982647.2655],
        [203867603.322, 163612745.223, 121765478.951, 84705239.0309, 54203841.1893, 34916553.3623],
        [193000258.246, 159633040.67, 120739509.738, 83739958.6425, 56950404.9026, 37220541.4696],
        [196964790.184, 156444685.295, 119223858.616, 85247685.164, 60502251.1293, 40693294.6208],
        [185401043.818, 147695755.167, 113467540.659, 87920455.11, 60709980.127, 42093478.8237],
        [167909073.769, 146713180.85, 116465787.385, 89571584.9871, 64926699.746, 43691869.1895],
        [167158116.581, 145482676.804, 114951518.904, 88475523.2607, 65030701.9542, 46120088.6093],
        [157760706.291, 137836455.92, 113379915.298, 91641365.1536, 68924700.2428, 48698033.483],
        [143607150.179, 127875584.892, 112520771.519, 93223216.4379, 70387887.7272, 50546712.6645],
        [112483697.906, 117150813.299, 107701424.781, 97751921.982, 81366149.7358, 64815248.0104],
        [98061096.2279, 107455892.672, 104297905.683, 100185926.349, 82831840.0415, 70943170.025],
        [89137895.3342, 102529524.376, 103609158.886, 100258967.138, 82953357.8248, 74711580.7698],
        [80061011.5245, 99701974.4374, 101135778.326, 103509334.667, 93903836.9186, 82575604.803],
        [76955476.7623, 99204914.2963, 100657431.412, 105555739.515, 97037484.1735, 86439016.9643],
        [69605396.4523, 94038539.7022, 99107027.4359, 106465379.024, 99860594.6244, 91198028.5089],
        [65217286.435, 91077737.2129, 97711464.8466, 106493430.55, 102497916.497, 94414439.934],
        [56404760.6472, 89515576.0369, 96725175.6269, 110780475.382, 106003429.523, 99779805.5687],
        [57623765.1379, 86000000.115, 92645901.5832, 109963760.39, 112042014.351, 105580958.636],
        [43952333.816, 79252453.6185, 86299264.174, 117994751.33, 118933764.105, 117810670.231],
        [24634373.903, 66065001.64, 68243438.02, 135426575.793, 137419433.105, 172318869.066],
        [18589595.0428, 57728776.2495, 62320493.2527, 139634562.057, 152773746.081, 223807728.344],
        [20594648.2386, 59271804.7584, 63480011.0556, 137984419.731, 170575636.626, 249672094.732],
        [19982219.0805, 51425904.8348, 62459589.4291, 133533308.408, 173764569.114, 340432929.048],
        [14233265.2396, 48314855.5797, 69876452.1042, 118939345.179, 202726701.441, 782283022.416],
        [12050167.2655, 41971380.0582, 67231664.3917, 117242812.869, 204972074.627, 1193155462.46],
        [217614427.542, 172388819.341, 121628347.075, 83816254.3477, 55033857.4819, 34150040.9715],
        [188956032.614, 152228600.204, 119348410.647, 87194122.6816, 60562480.3003, 39268083.3223],
        [181744158.147, 148070781.167, 115670547.145, 89247318.1735, 64040594.4522, 43685819.4112],
        [151073623.188, 144854097.547, 117506093.143, 91800688.1804, 68393724.3162, 47858125.7734],
        [104973486.386, 113192497.975, 107358293.398, 99298184.4158, 82854153.6826, 66448368.5807],
        [96800968.5366, 108813789.1, 107709745.619, 100869902.616, 87949796.8523, 73853744.7042],
        [77542089.7039, 100444892.881, 104405576.351, 103448058.566, 92715073.403, 80484729.3368],
        [73902932.3516, 98357005.817, 101596272.347, 106187268.512, 94949640.7439, 84977461.0913],
        [65723772.4606, 92220226.5048, 100949356.3, 107088570.17, 98364839.1687, 88907820.7861],
        [58227272.214, 88096043.6116, 97778537.1547, 109673326.262, 103072681.201, 95546867.7228],
        [52242736.0392, 79021365.3116, 92296514.2048, 110699401.452, 111373663.316, 104556413.811],
        [42387748.6585, 75036578.343, 88900875.6779, 113953224.692, 119871722.227, 118650222.998],
        [42552145.4285, 73277455.5577, 86638982.1563, 117509435.252, 122443225.523, 127794645.949],
        [29487207.2902, 62587561.9928, 80747503.8074, 122638642.662, 133371495.651, 146929399.381],
        [232212730.923, 172748200.813, 123164310.855, 83234536.1767, 52534357.4283, 30163940.3295],
        [209046363.073, 165889951.977, 122621699.858, 83335081.7221, 54798420.2942, 34266296.8977],
        [175059723.489, 158635872.558, 119903390.149, 84477402.7274, 59083638.5898, 39711091.9335],
        [180976714.043, 154287876.699, 115600774.565, 88075679.1989, 61996037.5575, 41220965.6938],
        [166541063.013, 144991639.006, 116876672.36, 86937102.9577, 64578343.8305, 44125027.459],
        [134791895.785, 127292591.696, 113389093.855, 94065035.7771, 71569091.3325, 50277016.4088],
        [116307285.53, 120873545.356, 110371130.191, 95960676.9701, 76880992.7215, 61613992.3649],
        [93643102.1393, 106738439.582, 105458128.829, 100696412.304, 87801584.4699, 73494131.1528],
        [70155926.5125, 92617201.6151, 99897419.2984, 107979933.81, 102133656.856, 91712156.019],
        [48432684.2664, 77640488.3084, 95885230.2566, 113033574.748, 113968272.876, 105856099.732],
        [44330245.1179, 76228877.0172, 92090755.5776, 115109943.263, 116827090.805, 111978280.037],
        [41986963.3733, 76941616.7393, 91077765.8441, 114116953.114, 120770471.227, 118076659.132],
        [24350677.148, 62062462.4579, 81391127.5714, 120199164.039, 139426704.127, 139859946.552],
        [16566555.6961, 51508541.8403, 74545214.652, 131282172.01, 163307573.327, 213211499.471],
        [15667048.0325, 54516055.2198, 85268892.8266, 137304083.086, 180118662.136, 244286666.947],
        [26899712.5025, 47612160.0504, 93480391.8501, 116366130.415, 245905247.734, 859551349.959],
        [77724573.3068, 24522031.096, 59371312.1115, 85723599.4897, 636830700.998, 5457044198.07],
        [184221616.692, 157528964.72, 121670502.315, 82992698.6482, 55715276.2256, 34907386.7376],
        [186240877.897, 146467471.52, 118213111.658, 86607956.5869, 59598204.4652, 37143728.2774],
        [136253992.631, 126603539.187, 114298371.853, 91497814.6939, 72230713.8519, 53987276.2465],
        [106314561.389, 113964562.672, 107918869.145, 100428267.401, 80796162.5065, 63503161.2242],
        [47356124.0695, 77562337.7631, 96053464.4237, 113941864.156, 112215268.984, 102384975.424],
        [23750150.5736, 55415541.0036, 87983441.2884, 123408442.537, 133986521.045, 127934773.286],
        [17255606.3062, 50453390.7123, 85294570.3812, 132676718.474, 154723573.633, 163970266.721],
        [11155571.6465, 38800946.9791, 78247709.2722, 140346162.545, 224168766.021, 484901098.72],
        [200559738.269, 159078279.738, 122953889.782, 82731894.4139, 55184352.7519, 34838689.9769],
        [165171567.434, 140314635.161, 117112512.952, 86173443.7999, 63763344.3168, 43066897.3621],
        [154588410.824, 135316209.281, 116328188.041, 86979136.0898, 64298781.8186, 45675697.948],
        [162002575.451, 135860600.864, 116417372.501, 87962745.6728, 67385382.8884, 47394394.2437],
        [128859771.183, 124587199.774, 112922590.213, 94837327.0756, 77299742.6425, 59873845.1031],
        [114276325.277, 122425104.132, 113704531.172, 95586394.9436, 81568052.7778, 66492550.8975],
        [69722356.9925, 92857115.5751, 107702036.034, 104686022.445, 90162037.5295, 75155699.3395],
        [53531426.5745, 86250467.0279, 100677800.573, 111415319.959, 103181772.44, 86715237.9212],
        [35730667.9461, 66247934.2, 96896061.1902, 120020733.258, 115654081.28, 99121977.5028],
        [26298204.4428, 62679511.988, 88463462.5452, 124553583.264, 136807712.722, 111284163.958],
        [15272331.1628, 44586299.5689, 78896976.0433, 131072821.633, 154407202.907, 143787641.777],
        [10706035.275, 38633176.148, 78176970.26, 143005159.177, 184412033.8, 235682483.354],
        [9664335.59595, 36576257.9713, 75767089.7455, 138404693.488, 226474614.578, 437137873.508],
        [83175694.5772, 88438166.4869, 94345688.9479, 103217951.644, 108653050.69, 31193044.2219],
        [60303529.1302, 74640785.895, 85751272.669, 113775336.998, 126310448.118, 143318527.736],
        [127382843.883, 127391061.218, 110617916.138, 84726610.041, 75952409.0468, 75952409.0468],
        [56804333.6506, 73257905.0076, 75274303.5987, 103673126.61, 105849015.291, 113943879.481],
        [36338997.8303, 65471934.6664, 80586441.3589, 124557270.448, 133214955.509, 184503981.16],
    ]

    return dict(zip(spectypes, flux_array))


def get_object_flux(
    spectype=None,
    teff=None,
    lumclass="V",
):

    spec_types = [
        "O5",
        "O9",
        "B0",
        "B1",
        "B3",
        "B6",
        "B8",
        "A0",
        "A2",
        "A5",
        "F0",
        "F2",
        "F5",
        "F8",
        "G0",
        "G2",
        "G5",
        "G8",
        "K0",
        "K2",
        "K5",
        "K7",
        "M0",
        "M2",
        "M4",
        "M5",
    ]

    if spectype is None:

        if teff is None:
            print("MUST PROVIDE SPECTRAL TYPE OR TEMPERATURE.\n\n")
        if np.isnan(teff):
            teff = 5800.0

        teff_ranges = np.array(
            [
                55_000,
                38_000,
                30_000,
                23_000,
                18_000,
                15_000,
                12_000,
                10_000,
                9_000,
                8300,
                7400,
                7000,
                6700,
                6300,
                6000,
                5800,
                5660,
                5440,
                5240,
                4960,
                4400,
                4000,
                3750,
                3600,
                3400,
                3200,
            ]
        )

        spectype = spec_types[np.argmin(np.abs(teff - teff_ranges))]
        flux_key = spectype + lumclass

    elif any(np.isin([spectype], ["Seyfert1", "Seyfert2", "QSO", "NGC1068", "Liner"])):
        flux_key = spectype

    else:
        flux_key = spectype

    flux_dict = get_object_flux_dictionary()

    if not (any(np.isin([flux_key], list(flux_dict.keys())))):
        print(flux_key + " NOT IN LIST OF SPECTRAL TYPES. PLEASE CHOOSE FROM BELOW: ")
        print(flux_dict.keys())

    return flux_dict[flux_key]


def calc_n_photon_per_sec(
    Vmag,
    spectype=None,
    teff=None,
    lumclass="V",
    fiber_setup="100",
    airmass=1.5,
    seeing=1.0,
):

    if fiber_setup == "300":
        # arguments for Low-Res mode.
        pinhole = 2.3
        binning = 5.0
        npix = 35.0 * 3.0
        m1size = 8.4  # // effic 2019
        insttrans = np.array([1.419, 2.064, 4.454, 3.982, 6.469, 8.261])

    elif fiber_setup == "200":
        pinhole = 1.5
        binning = 2.0
        npix = 17.0 * 5.0
        m1size = 8.4  # // effic 2019
        # insttrans=[0.90,2.00,3.83,4.73,4.93,6.07];
        insttrans = [0.545, 1.697, 2.528, 3.554, 4.065, 4.045]

    elif fiber_setup == "100":
        pinhole = 0.74
        binning = 1.0
        npix = 9.0 * 7.0
        m1size = 8.4  # effic 2019
        # insttrans=[0.88,1.35,2.40,3.98,2.84,4.28];
        # effic 2021
        insttrans = [0.448, 2.015, 1.992, 1.637, 2.278, 0.986]

    elif fiber_setup == "VATT":
        pinhole = 7.5
        binning = 1.0
        npix = 20.0 * 9.0
        m1size = 1.9
        insttrans = [0.00, 0.50, 1.50, 4.00, 5.00, 4.00]
        binocular = False

    elif fiber_setup == "POL":
        pinhole = 1.5
        binning = 2.0
        npix = 17.0 * 5.0
        m1size = 8.4
        insttrans = [0.39, 1.34, 5.10, 3.87, 7.13, 7.61]

    else:
        print("fiber_setup must be '300', '200', '100', 'VATT' or 'POL' ")

    # properties for each cross-disperser
    teltrans_arr = 0.01 * np.array(insttrans)
    flam_arr = [0.402704, 0.450081, 0.505877, 0.582962, 0.687765, 0.827176]
    fwhm_arr = np.array([0.000698, 0.000781, 0.000887, 0.000990, 0.001182, 0.001459])
    readnoise_arr = [3.5, 3.5, 3.5, 4.0, 4.0, 4.0]

    # Add extinction to Vmag, calculate incoming flux
    kparr = np.array([0.3035, 0.2067, 0.1487, 0.1180, 0.0678, 0.0418])
    vmag_extincted = Vmag + kparr * airmass  # extinction through atmosphere
    zeropoint_fluxarr = get_object_flux(spectype=spectype, teff=teff, lumclass=lumclass)

    # transmission through fiber pinhole
    trans = np.power(math.erf(pinhole / (seeing * 1.202)), 2)
    area = np.power(m1size, 2) * np.pi / 4.0

    n_photons_per_sec = (
        zeropoint_fluxarr * fwhm_arr * binning * area * teltrans_arr * trans * np.power(10, -0.4 * vmag_extincted)
    )

    return n_photons_per_sec, npix


def pepsi_snr(
    Vmag, exptime, spectype=None, teff=None, lumclass="V", fiber_setup="100", airmass=1.5, seeing=1.0, binocular=False
):
    """
    Vmag: V-band magnitude (Vega) of your source
    exptime: Exposure time in seconds
    spectype: Spectral type of stellar source, also can be any of ['Seyfert1', 'Seyfert2', 'QSO', 'NGC1068', 'Liner']
    teff: effective temperature in Kelvin, can provide instead of Spectral Type
    resolution: 'low' (R=50,000), 'med' (R=130,000), or 'high' (R=250,000)
    binocular: True/False
    airmass: observed airmass
    seeing: observed seeing (FWHM in arcsec)
    """

    if fiber_setup == "VATT":
        binocular = False

    n_photon_per_sec, npix = calc_n_photon_per_sec(Vmag, spectype, teff, lumclass, fiber_setup, airmass, seeing)

    # Readnoise for the Blue and Red CCD, respectively
    readnoise_arr = [3.5, 3.5, 3.5, 4.0, 4.0, 4.0]

    n_photon = n_photon_per_sec * exptime
    snratio = n_photon / np.sqrt(n_photon + npix * np.power(readnoise_arr, 2.0))

    if binocular:
        snratio *= np.sqrt(2.0)

    return np.round(snratio, 1)


def pepsi_exptime(
    Vmag, snratio, spectype=None, teff=None, lumclass="V", fiber_setup="100", airmass=1.5, seeing=1.0, binocular=False
):
    """
    snr: desired signal-to-noise ratio
    Vmag: V-band magnitude (Vega) of your source
    spectype: Spectral type of stellar source, also can be any of ['Seyfert1', 'Seyfert2', 'QSO', 'NGC1068', 'Liner']
    teff: effective temperature in Kelvin, can provide instead of Spectral Type
    resolution: 'low' (R=50,000), 'med' (R=130,000), or 'high' (R=250,000)
    binocular: True/False
    air: observed airmass
    see: observed seeing (FWHM in arcsec)
    """

    n_photon_per_sec, npix = calc_n_photon_per_sec(Vmag, spectype, teff, lumclass, fiber_setup, airmass, seeing)

    if fiber_setup == "VATT":
        binocular = False

    if binocular:
        snratio /= np.sqrt(2)

    # Readnoise for the Blue and Red CCD, respectively
    readnoise_arr = [3.5, 3.5, 3.5, 4.0, 4.0, 4.0]

    # Solve array equation for texp:
    # S/N = n_photon_per sec * texp / sqrt( n_photon_per_sec * texp + npix*readnoise^2 )
    a = -np.power(n_photon_per_sec, 2)
    b = np.power(snratio, 2) * n_photon_per_sec
    c = np.power(snratio, 2) * npix * np.power(readnoise_arr, 2.0)

    exptime = (-b - np.sqrt(b**2.0 - 4 * a * c)) / (2 * a)

    return np.round(exptime, 1)


def add_pepsi_params(
    df: pd.DataFrame, fiber: str, cd_blue: int, cd_red: int, snr: int, binocular: bool = False, **kwargs
) -> pd.DataFrame:
    df["pepsi_fiber"] = fiber
    df["pepsi_cd_blue"] = cd_blue
    df["pepsi_cd_blue_num_exp"] = 1
    df["pepsi_cd_red"] = cd_red
    df["pepsi_cd_red_num_exp"] = 1
    df["pepsi_snr"] = snr
    df["pepsi_exp_time"] = [
        pepsi_exptime(vmag, snr, teff=teff, fiber_setup=fiber, binocular=binocular)[cd_red - 1]
        for vmag, teff in df[["Vmag", "Teff"]].values
    ]
    df["pepsi_priority"] = ""
    df["pepsi_notes"] = ""
    return df


from astropy.coordinates import Angle
from astropy.table import Table
import astropy.units as u
import numpy as np
import pandas as pd


def write_lbt_readme_file(file_base: str, targets: pd.DataFrame) -> str:
    # build up a new table in proper format
    readme = Table()
    readme["Target Name"] = targets["target_name"]
    readme["ra"] = targets["ra"]
    readme["RA"] = [Angle(ra, unit=u.deg).to_string(unit=u.hour, decimal=False, precision=2, sep=":") for ra in targets["ra"]]
    readme["DEC"] = [
        Angle(dec, unit=u.deg).to_string(unit=u.deg, decimal=False, precision=2, sep=":", alwayssign=True)
        for dec in targets["dec"]
    ]
    readme["Vmag"] = targets["Vmag"]
    readme["Teff"] = targets["Teff"]
    readme["Fiber"] = targets["pepsi_fiber"]
    readme["BLUE Cross Disperser"] = targets["pepsi_cd_blue"]
    readme["BLUE CD NExp"] = targets["pepsi_cd_blue_num_exp"]
    readme["BLUE CD Exp Time"] = targets["pepsi_exp_time"]
    readme["RED Cross Disperser"] = targets["pepsi_cd_red"]
    readme["RED CD NExp"] = targets["pepsi_cd_red_num_exp"]
    readme["RED CD Exp Time"] = targets["pepsi_exp_time"]
    readme["Desired BLUE SNR"] = targets["pepsi_snr"]
    readme["Desired RED SNR"] = targets["pepsi_snr"]
    readme["Priority"] = targets["pepsi_priority"]
    readme["Notes"] = targets["pepsi_notes"]
    readme.sort("ra")
    readme.remove_column("ra")

    # save target list as csv
    readme.write(file_base + ".csv")

    # make the readme file by prepending/appending the header/footer info
    readme_header = open(file_base + ".README.header", "r").readlines()
    readme_footer = open(file_base + ".README.footer", "r").readlines()
    output = ""
    for line in readme_header + readme.pformat_all() + readme_footer:
        output += line.rstrip() + "\n"
    with open(file_base + ".README", "w") as f:
        f.write(output)
    return output
