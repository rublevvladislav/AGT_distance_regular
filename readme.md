# Дистранционно регулярные графы
### Требования
Требования указаны в файле requirements.txt
для установки зависимостей можно использовать команду: 
> pip install -r requirements.txt
### Основные функции
** Являтеся ли граф дистанционно регулярным? ** - функция is_distance_regular_graph(adjacency_matrix) 
> принимает на вход матрицу смежности
> возвращает False или [True, массив пересечений] 
** Нахождение массива пересечений ** - функция intersect_array(adjacency_matrix)
> принимает на вход матрицу смежности
> возвращает массив пересечений
** Нахождение слоевого представления ** - функция layed_represent(start_point, adjacency_matrix)
> принимает на вход матрицу смежности и начальную точку для отсчета расстояний
> возвращает слоевое представление в виде массива
** Вычисление спектра графа (собственные числа и их кратности) ** - функция get_spectrum(adjacency_matrix)
> принимает на вход матрицу смежности
> возвращает два массива: собственных значений и кратностей