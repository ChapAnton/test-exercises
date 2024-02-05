from orbital_params import construct_orbital_params
import numpy as np
from physical_const import *


# начальный момент времени в юлианских днях
JD_START = 8084.05644194318
# количество секунд в сутках
SECONDS_PER_DAY = 86400


def count_position(semi_major_axis, ecc_anomaly, eccentricity, transf_mat):
    """
    расчет положения спутника
    Returns:
        pos (dict): 
            x (float): x - компонента радиус-вектора спутника в осях ИСК
            y (float): y - компонента радиус-вектора спутника в осях ИСК
            z (float): z - компонента радиус-вектора спутника в осях ИСК
    """
    x, y, z = (semi_major_axis * (transf_mat[0]
                                  * (np.cos(ecc_anomaly) - eccentricity) + transf_mat[1] *
                                  np.sqrt(1 - eccentricity ** 2) * np.sin(ecc_anomaly)))
    pos = {'x': float(x), 'y': float(y), 'z': float(z)}

    return pos


def count_velocity(semi_major_axis, ecc_anomaly, eccentricity, transf_mat):
    """
    расчет вектора скорости спутника 
    Returns:
        vel (dict): 
            vel_x (float): x - компонента скорости спутника в осях ИСК
            vel_y (float): y - компонента скорости спутника в осях ИСК
            vel_z (float): z - компонента скорости спутника в осях ИСК
    """
    vel_x, vel_y, vel_z = (np.sqrt(MU / semi_major_axis)
                           / (1 - eccentricity * np.cos(ecc_anomaly))) * ((- transf_mat[0]) *
                                                                          np.sin(ecc_anomaly) + transf_mat[1] * (np.sqrt(1 - eccentricity ** 2)
                                                                                                                 * np.cos(ecc_anomaly)))
    vel = {'vel_x': float(vel_x), 'vel_y': float(vel_y), 'vel_z': float(vel_z)}

    return vel


class Satellite():
    """
    класс описывает параметры движения спутника
    """

    def __init__(self, start_pos, start_vel, time_begin, target_point):
        """
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
            time_begin (float): начальный момент времени
            target_point (dict) : координаты целевой точки наведения в ИСК, 
            вида {'x':x, 'y':y, 'z': z}, где
                x (float): x - координата целевой точки
                y (float): y - координата целевой точки
                z (float): z - координата целевой точки
        """
        self._orb_params = construct_orbital_params(start_pos, start_vel)
        self._transf_mat = None
        self._time_begin = time_begin
        self._current_time = time_begin
        self.pos = {'x': None, 'y': None, 'z': None}
        self.vel = {'vel_x': None, 'vel_y': None, 'vel_z': None}
        self._target_point = target_point

    def update_params(self, time):
        """
        расчет векторов положения и скорости спутника на момент времени time
        Args:
            time (float): время в секундах  
        """
        self._orb_params.update_anomaly(time)
        self._current_time = time
        self.pos = count_position(self._orb_params._semi_major_axis,
                                  self._orb_params._ecc_anomaly, self._orb_params._eccentricity, self._transf_mat)

        self.vel = count_velocity(self._orb_params._semi_major_axis,
                                  self._orb_params._ecc_anomaly, self._orb_params._eccentricity, self._transf_mat)

    def get_params(self):
        """
        возврат вектора состояния спутника
        Returns:
            result (list): вектор состояния спутника в момент времени time вида
            [julian_date, x, y, z, vel_x, vel_y, vel_z], где
            julian_date (float): время в формате количества юлинский дней в шкале TT, прошедших от эпохи J2000 
            x (float): x - компонента радиус-вектора спутника в осях ИСК
            y (float): y - компонента радиус-вектора спутника в осях ИСК
            z (float): z - компонента радиус-вектора спутника в осях ИСК
            vel_x (float): x - компонента скорости спутника в осях ИСК
            vel_y (float): y - компонента скорости спутника в осях ИСК
            vel_z (float): z - компонента скорости спутника в осях ИСК
        """
        julian_date = JD_START + \
            (self._time_begin + self._current_time) / SECONDS_PER_DAY
        result = [julian_date, self.pos['x'], self.pos['y'], self.pos['z'],
                  self.vel['vel_x'], self.vel['vel_y'], self.vel['vel_z']]

        return result


def construct_satellite(start_pos, start_vel, time_begin, target_point):
    """
    создание и инициализация объекта Satellite
    Args:
       start_pos (dict): словарь, содержащий начальное положение спутника
            в ИСК, вида {'x': x, 'y': y,'z': z}, где
                x (float): x - компонента радиус-вектора спутника 
                y (float): y - компонента радиус-вектора спутника 
                z (float): z - компонента радиус-вектора спутника 
        start_velocity (list): словарь, содержащий начальный вектор скорости 
        спутника в ИСК, вида {'vel_x': vel_x, 'vel_y': vel_y,'vel_z': vel_z}, 
        где 
            vel_x (float): x - компонента скорости спутника
            vel_y (float): y - компонента скорости спутника
            vel_z (float): z - компонента скорости спутника
        time_begin (float): начальный момент времени
    Returns:
        satellite(obj): экземляр класса Satellite
    """
    satellite = Satellite(start_pos, start_vel, time_begin, target_point)
    ascending_node = satellite._orb_params._ascending_node
    periapsis_arg = satellite._orb_params._periapsis_arg
    inclination = satellite._orb_params._inclination

    m11 = np.cos(ascending_node) * np.cos(periapsis_arg) - \
        np.sin(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)
    m21 = np.sin(ascending_node) * np.cos(periapsis_arg) + \
        np.cos(ascending_node) * np.sin(periapsis_arg) * np.cos(inclination)
    m31 = np.sin(periapsis_arg) * np.sin(inclination)

    m12 = - np.cos(ascending_node) * np.sin(periapsis_arg) - \
        np.sin(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)
    m22 = - np.sin(ascending_node) * np.sin(periapsis_arg) + \
        np.cos(ascending_node) * np.cos(periapsis_arg) * np.cos(inclination)
    m32 = np.cos(periapsis_arg) * np.sin(inclination)

    satellite._transf_mat = np.array([[m11, m21, m31], [m12, m22, m32]])

    return satellite
