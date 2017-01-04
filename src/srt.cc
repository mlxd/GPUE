
double sepAvg(int2 *vArray, int2 centre, int length){
	double result=0.0;// = sqrt( pow(centre.x - v_array[0].x,2) + pow(centre.y - v_array[0].y,2));
	for (int i=0; i<length; ++i){
		result += sqrt( pow(centre.x - v_array[i].x,2) + pow(centre.y - v_array[i].y,2));
	}
	return result/length;
}
