import math
import numpy as np
import time

# Функции преобразования координат

def polar_to_cartesian(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def cartesian_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def spherical_to_cartesian(r, theta, phi):
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z

def cartesian_to_spherical(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y, x)
    phi = np.arccos(z / r)
    return r, theta, phi

def distance_cartesian_2d(p1, p2):
    return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

def distance_polar(p1, p2):
    r1, theta1 = p1
    r2, theta2 = p2
    return np.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * np.cos(theta2 - theta1))

def distance_cartesian_3d(p1, p2):
    return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)

def distance_spherical(p1, p2):
    r1, theta1, phi1 = p1
    r2, theta2, phi2 = p2
    return np.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * (np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2) + np.cos(theta1) * np.cos(theta2)))

def great_circle_distance(p1, p2, radius):
    theta1, phi1 = p1[1], p1[2]
    theta2, phi2 = p2[1], p2[2]
    return radius * np.arccos(np.sin(theta1) * np.sin(theta2) + np.cos(theta1) * np.cos(theta2) * np.cos(phi2 - phi1))

def main():
    # Генерация точек
    n_points = 10000
    r_values = np.random.default_rng().uniform(0.0, 10.0, n_points)
    theta_values = np.random.default_rng().uniform(0, np.pi, n_points)
    phi_values = np.random.default_rng().uniform(0, np.pi, n_points)

    polar_coords = np.column_stack((r_values, theta_values))
    spherical_coords = np.column_stack((r_values, theta_values, phi_values))

    # 2D: Полярные -> Декартовы -> Полярные
    cartesian_coords = np.array([polar_to_cartesian(r, theta) for r, theta in polar_coords])
    polar_coords_back = np.array([cartesian_to_polar(x, y) for x, y in cartesian_coords])

    # 3D: Сферические -> Декартовы -> Сферические
    cartesian_coords_3d = np.array([spherical_to_cartesian(r, theta, phi) for r, theta, phi in spherical_coords])
    spherical_coords_back = np.array([cartesian_to_spherical(x, y, z) for x, y, z in cartesian_coords_3d])

    # Проверка корректности
    assert np.allclose(polar_coords, polar_coords_back), "Ошибка в преобразовании координат (2D)"
    assert np.allclose(spherical_coords, spherical_coords_back), "Ошибка в преобразовании координат (3D)"

    # Бенчмаркинг
    start_time = time.time()
    for i in range(n_points - 1):
        d_cartesian_2d = distance_cartesian_2d(cartesian_coords[i], cartesian_coords[i + 1])
    end_time = time.time()
    print(f"Время выполнения (Декартовые 2D): {end_time - start_time:.4f} секунд")

    start_time = time.time()
    for i in range(n_points - 1):
        d_cartesian_3d = distance_cartesian_3d(cartesian_coords_3d[i], cartesian_coords_3d[i + 1])
    end_time = time.time()
    print(f"Время выполнения (Декартовые 3D): {end_time - start_time:.4f} секунд")

    start_time = time.time()
    for i in range(n_points - 1):
        d_polar = distance_polar(polar_coords[i], polar_coords[i + 1])
    end_time = time.time()
    print(f"Время выполнения (Полярные): {end_time - start_time:.4f} секунд")

    start_time = time.time()
    for i in range(n_points - 1):
        d_spherical = distance_spherical(spherical_coords[i], spherical_coords[i + 1])
    end_time = time.time()
    print(f"Время выполнения (Объем сферы): {end_time - start_time:.4f} секунд")

    start_time = time.time()
    for i in range(n_points - 1):
        d_great_circle = great_circle_distance(spherical_coords[i], spherical_coords[i + 1], r_values[i])
    end_time = time.time()
    print(f"Время выполнения (Плоскость сферы): {end_time - start_time:.4f} секунд")

if __name__ == "__main__":
    main()