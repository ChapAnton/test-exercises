from satellite import construct_satellite
from pymap3d import ecef2eci, geodetic2ecef
import datetime


# начальный момент времени в секундах
TIME_BEGIN = 0 
# конечный момент времени в секундах
TIME_END = 20000
# временной шаг
STEP = 1
# начальное положение спутника в ИСК
START_POS = {'x': 4362521.19692133, 'y': -2174459.71448059, 
            'z': 4720847.40402189}
# начальный вектор скорости спутника в ИСК          
START_VEL = {'vel_x': 5356.39915538069,
            'vel_y': 4741.41348686709, 'vel_z': -2761.5472632395}
# координаты целевой точки на поверхности Земли
TARGET_POINT = {'lat': 45.920266, 'lon': 63.342286}
# начальный момент времени в календарном формате
CALENDAR_TIME_START = datetime.datetime(2022, 2, 18, 13, 21, 16, 584) 


def convert_geo2eci(target_point, time):
    """
    преобразование координат целевой точки из геодезической СК в ИСК
    Args:
        target_point (dict): координатцы целевой точки в геодезической СК
        time (datetime): время наблюдения

    Returns:
       target_point_eci (dict): словарь, содержащий координаты целевой точки 
       в ИСК вида {'x': x, 'y': y, 'z': z}, где
            x (float): x - компонента радиус-вектора целевой точки 
            y (float): y - компонента радиус-вектора целевой точки 
            z (float): z - компонента радиус-вектора целевой точки 
    """
    x_ecef, y_ecef, z_ecef = geodetic2ecef(target_point['lat'], target_point['lon'], 0)
    x_eci, y_eci, z_eci = ecef2eci([x_ecef, y_ecef, z_ecef], time)
    target_point_eci = {'x': x_eci, 'y': y_eci, 'z': z_eci}
    return target_point_eci


def process(time_begin, time_end, step=1):
    """
    запись в файл вектора состояния спутника в заданном промежутке времени
    Args:
        time_begin (float): начальный момент времени
        time_end (_type_): конечный момент времени
        step (int, optional): временной шаг
    """
    target_point = convert_geo2eci(TARGET_POINT, CALENDAR_TIME_START)
    satellite = construct_satellite(START_POS, START_VEL, TIME_BEGIN, target_point)
    with open('result.txt', 'w') as f:
        for i in range(time_begin, time_end, step):
            satellite.update_params(i)
            out_data = satellite.get_params()
            out_data = ' '.join(str(i) for i in out_data) + '\n'
            f.writelines(out_data)
            
if __name__ == "__main__":
    process(TIME_BEGIN, TIME_END, STEP)
