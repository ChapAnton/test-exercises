import numpy as np
import kepler
from physical_const import *


class OrbitalParams():
    """
    класс описывает кеплеровы параметры орбиты
    """
    def __init__(self):
        self._semi_major_axis = None
        self._eccentricity = None
        self._inclination = None
        self._ascending_node = None
        self._true_anomaly = None
        self._periapsis_arg = None
        self._mean_anomaly_begin = None
        self._average_distance = None
        self._ecc_anomaly = None
        self._mean_anomaly = None

    def update_anomaly(self, time):
        """
        Обновление аномалий орбиты по времени 

        Args:
            time (float): момент времени time
        """
        self._mean_anomaly = self._mean_anomaly_begin + time \
            * self._average_distance
        self._ecc_anomaly, cos_true_anomaly, _ = kepler.kepler(
            self._mean_anomaly, self._eccentricity)
        self._true_anomaly = np.arccos(cos_true_anomaly)


def construct_orbital_params(start_pos, start_vel):
    """
    создание и инициализация объекта OrbitalParams

    Args:
        start_pos (dict): словарь, содержащий начальное положение спутника
        в ИСК, вида {'x': x, 'y': y,'z': z}, где
            x (float): x - компонента радиус-вектора спутника 
            y (float): y - компонента радиус-вектора спутника 
            z (float): z - компонента радиус-вектора спутника 
        start_vel (dict): словарь, содержащий начальный вектор \
        скорости спутника в ИСК, вида {'vel_x': vel_x, 'vel_y': vel_y,
        'vel_z': vel_z}, где 
            vel_x (float): x - компонента скорости спутника 
            vel_y (float): y - компонента скорости спутника 
            vel_z (float): z - компонента скорости спутника 
    Returns:
        orb_params(obj): экземпляр класса OrbitalParams
    """
    orb_params = OrbitalParams()

    pos_proj = list(start_pos.values())
    vel_proj = list(start_vel.values())

    pos_norm = np.linalg.norm(pos_proj, ord=2)
    vel_norm = np.linalg.norm(vel_proj, ord=2)

    orb_moment = count_orb_moment(pos_proj, vel_proj)
    orb_moment_norm = np.linalg.norm(orb_moment)

    keplerian_energy = count_keplerian_energy(pos_norm, vel_norm)

    focal_distance = count_focal_distance(orb_moment_norm)

    orb_params._semi_major_axis = count_semi_major_axis(keplerian_energy)

    orb_params._eccentricity = count_eccentricity(
        keplerian_energy, orb_moment_norm)

    orb_params._true_anomaly = count_true_anomaly(
        pos_proj, vel_proj, orb_moment_norm, pos_norm, focal_distance)

    orb_params._periapsis_arg = count_periapsis_arg(
        start_pos, orb_moment, orb_moment_norm, orb_params._true_anomaly)

    orb_params._inclination = count_inclination(
        orb_moment, orb_moment_norm)

    orb_params._ascending_node = count_ascending_node(orb_moment)

    ecc_anomaly_begin = count_ecc_anomaly_begin(
        pos_proj, vel_proj, orb_params._semi_major_axis, pos_norm)

    orb_params._mean_anomaly_begin = count_mean_anomaly_begin(
        ecc_anomaly_begin, pos_proj, vel_proj, orb_params._semi_major_axis)

    orb_params._average_distance = count_average_distance(
        orb_params._semi_major_axis)

    orb_params._mean_anomaly = orb_params._mean_anomaly_begin
    orb_params._ecc_anomaly = ecc_anomaly_begin

    return orb_params


def count_keplerian_energy(pos_norm, vel_norm):
    """ 
    расчет энергии кеплерова движения
    Args:
        pos_norm (float): евклидова норма радиус вектора спутника 
            в начальный момент времени
        vel_norm (float): евклидова норма  вектора скорости спутника
            в начальный момент времени
    Returns:
        keplerian_energy (float): энергия кеплерова движения
    """
    keplerian_energy = ((vel_norm ** 2) / 2) - (MU / pos_norm)

    return keplerian_energy


def count_orb_moment(pos_proj, vel_proj):
    """
    расчет орбитального момента
    Args:
        pos_proj (list): проекции вектора положения спутника на оси ИСК
        vel_proj (list): проекции вектора скорости спутника на оси ИСК
    Returns:
        orb_moment (ndarray): орбитальный момент 
    """
    orb_moment = np.cross(pos_proj, vel_proj)

    return orb_moment


def count_eccentricity(keplerian_energy, orb_moment_norm):
    """
    расчет эксцентрисистета
    Args:
        keplerian_energy (float): энергия кеплерова движения
        orb_moment_norm (float): евклидова норма вектора орбитального момента 
    Returns:
        eccentricity (float): эксцентрисистет орбиты
    """
    eccentricity = np.sqrt(
        (1 + (2 * keplerian_energy) * (orb_moment_norm ** 2 / MU ** 2)))

    return eccentricity


def count_true_anomaly(pos_proj, vel_proj, orb_moment_norm, pos_norm, focal_distance):
    """
    расчет истинной аномалии
    Args:
        pos_proj (list): проекции вектора положения спутника на оси ИСК
        vel_proj (list): проекции вектора скорости спутника на оси ИСК
        orb_moment_norm (float): евклидова норма вектора орбитального момента 
        pos_norm (float): евклидова норма радиус вектора спутника
        focal_distance (float): фокальное расстояние
    Returns:
        true_anomaly (float): истинная аномалия
    """
    true_anomaly = np.arctan2(np.dot(pos_proj, vel_proj) 
        / orb_moment_norm, 1 - pos_norm / focal_distance)

    return true_anomaly


def count_periapsis_arg(start_pos, orb_moment, orb_moment_norm, true_anomaly):
    """
    расчет аргумента перицентра
    Args:
        start_pos (dict): словарь, содержащий начальное положение спутника
        в ИСК, вида {'x': x, 'y': y,'z': z}, где
            x (float): x - компонента радиус-вектора спутника 
            y (float): y - компонента радиус-вектора спутника 
            z (float): z - компонента радиус-вектора спутника 
        orb_moment (ndarray): орбитальный момент 
        orb_moment_norm (float): евклидова норма вектора орбитального момента 
        true_anomaly (float): истинная аномалия
    Returns:
        periapsis_arg (float): аргумент перицентра
    """
    periapsis_arg = np.arctan2(start_pos['z'], (start_pos['y'] * orb_moment[0] -
        start_pos['x'] * orb_moment[1]) / orb_moment_norm) - true_anomaly

    return periapsis_arg


def count_semi_major_axis(keplerian_energy):
    """
    расчет большой полуоси 
    Args:
        keplerian_energy (float): энергия кеплерова движения
    Returns:
        semi_major_axis (float): большая полуось 
    """
    semi_major_axis = - MU / (2 * keplerian_energy)

    return semi_major_axis


def count_focal_distance(orb_moment_norm):
    """
    расчет фокального расстояния
    Args:
        orb_moment_norm (float): евклидова норма вектора орбитального момента 
    Returns:
        focal_distance (float): фокальное расстояние
    """
    focal_distance = orb_moment_norm ** 2 / MU

    return focal_distance


def count_inclination(orb_moment, orb_moment_norm):
    """
    расчет наклонения орбиты
    Args:
        orb_moment (ndarray): орбитальный момент 
        orb_moment_norm (float): евклидова норма вектора орбитального момента
    Returns:
        inclination (float): наклонение орбиты 
    """
    inclination = np.arccos(orb_moment[2] / orb_moment_norm)

    return inclination


def count_ascending_node(orb_moment):
    """
    расчет долготы восходящего узла
    Args:
        orb_moment (ndarray): орбитальный момент 
    Returns:
        ascending_node (float): долгота восходящего узла
    """
    ascending_node = np.arctan2(orb_moment[0], - orb_moment[1])

    return ascending_node


def count_ecc_anomaly_begin(pos_proj, vel_proj, semi_major_axis, pos_norm):
    """
    расчет эксцентрической аномалии на начальный момент времени
    Args:
        pos_proj (list): проекции вектора положения спутника на оси ИСК
        vel_proj (list): проекции вектора скорости спутника на оси ИСК
        semi_major_axis (float): большая полуось 
        pos_norm (float): евклидова норма радиус вектора спутника
    Returns:
        ecc_anomaly_begin (float): эксцентрическая аномалия на начальный момент времени
    """
    count_ecc_anomaly_begin = np.arctan2(np.dot(pos_proj, vel_proj) / 
        (np.sqrt(MU * semi_major_axis)), 1 - pos_norm / semi_major_axis)

    return count_ecc_anomaly_begin


def count_mean_anomaly_begin(ecc_anomaly_begin, pos_proj, vel_proj, semi_major_axis):
    """
    расчет средней аномалии на начальный момент времени
    Args:
        ecc_anomaly_begin (float): эксцентрическая аномалия на начальный момент времени
        pos_proj (list): проекции вектора положения спутника на оси ИСК
        vel_proj (list): проекции вектора скорости спутника на оси ИСК
        semi_major_axis (float): большая полуось 

    Returns:
        mean_anomaly_begin (float): средняя аномалия на начальный момент времени
    """
    mean_anomaly_begin = ecc_anomaly_begin - \
        np.dot(pos_proj, vel_proj) / (np.sqrt(MU * semi_major_axis))

    return mean_anomaly_begin


def count_average_distance(semi_major_axis):
    """
    расчет среднего движения
    Args:
        semi_major_axis (float): большая полуось 

    Returns:
        average_distance (float): среднее движение
    """
    average_distance = np.sqrt(MU / semi_major_axis ** 3)

    return average_distance
