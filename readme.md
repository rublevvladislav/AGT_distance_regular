# Дистранционно регулярные графы
### Требования
Программа написана на языке python 3.6
Можно использовать как самостоятельный скрипт или подключать в виде библиотеки в другой проект.
Зависимости указаны в файле requirements.txt,
для их установки можно использовать команду: 
> pip install -r requirements.txt
### Класс Graph
Есть два варианта работы с классом Graph: 
1) Передать конструктору матрицу смежности
2) Оставить конструктор пустым и использовать функции add_vertex() и add_edge()
### Основные функции класса Graph
**Являтеся ли граф дистанционно регулярным?** - функция is_distance_regular_graph() 
- *принимает на вход матрицу смежности (либо передана конструктору, либо автоматически сгенерирована)*
- *возвращает False или [True, массив пересечений]* 

**Нахождение массива пересечений** - функция intersect_array()
- *принимает на вход матрицу смежности (либо передана конструктору, либо автоматически сгенерирована)*
- *возвращает массив пересечений*

**Нахождение слоевого представления** - функция layed_represent(start_point)
- *принимает на вход матрицу смежности и начальную точку для отсчета расстояний*
- *возвращает слоевое представление в виде массива*

  функция draw_layer_represent(start_point)
- *принимает на вход матрицу смежности и начальную точку для отсчета расстояний*
- *возвращает слоевое представление в виде графика*

**Вычисление спектра графа (собственные числа и их кратности)** - функция get_spectrum()
- *принимает на вход матрицу смежности*
- *возвращает два массива: собственных значений и кратностей*

### Ограничения
-Граф должен быть неориентированным