`timescale 1ns/10ps

module adc_dumper(ok_to_sample,
                         dig_l,
                         dig_r,
			 ana_clk, //rising edge when data is stable
			 ana_l,
			 ana_r);

input ok_to_sample;
input [19:0] dig_l;
input [19:0] dig_r;
input ana_clk;
input ana_l;
input ana_r;


integer file_l_id;
integer file_r_id;
integer file_la_id;
integer file_ra_id;


initial file_l_id = $fopen ("/local/scratch/simulation/adc_analog_2channel_simimage/left.txt");
initial file_r_id = $fopen ("/local/scratch/simulation/adc_analog_2channel_simimage/right.txt");
initial file_la_id = $fopen ("/local/scratch/simulation/adc_analog_2channel_simimage/left_a.txt");
initial file_ra_id = $fopen ("/local/scratch/simulation/adc_analog_2channel_simimage/right_a.txt");

always@(posedge ok_to_sample)
 $fdisplay(file_l_id,"%d",dig_l);

always@(posedge ok_to_sample)
 $fdisplay(file_r_id,"%d",dig_r);

always@(posedge ana_clk)
 $fdisplay(file_la_id,"%d",ana_l);

always@(posedge ana_clk)
 $fdisplay(file_ra_id,"%d",ana_r);

endmodule 
