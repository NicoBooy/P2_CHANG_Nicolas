# P2_CHANG_Nicolas

###### In that context, we will only consider constraints with "<=" and x,y >= 0

The instruction structure is like this :

For example : 
<pre markdown="1">  --obj_x 1 --obj_y 2 --maximize ` 
  --constraint '-1' '1' '2' ` 
  --constraint '1' '2' '8' ` 
  --constraint '1' '0' '6'  </pre>

is the traduction of : 
```
Objective function : z = x + 2y
subject to :
- -x + y <= 2
- x + 2y <= 8
- x <= 6
```
